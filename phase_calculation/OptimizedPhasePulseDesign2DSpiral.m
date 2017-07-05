% Script for 2D RF pulse design for calculating optimized target phase with
% certain acceleration factors for VCC reconstruction, as used in the
% manuscript (vendor-specific data loading and waveform formatting is
% omitted)
%
% Object phase is calculated after coil combination, latter done with ESPIRIT.
% Pulse design is done by using the method described in:
% Iterative RF Pulse Design for Multidimensional, Small-Tip-Angle Selective Excitation
% Chun-yu Yip, Jeffrey A. Fessler, and Douglas C. Noll, MRM 2005
%
% now it uses spiral excitation kspace traverse
%
% Written by Adam Kettinger, Budapest, Hungary, 2016. 03. 16. and later on
%
% see corresponding manuscript in References of the Readme.
%
% some housekeeping and initialization:
close all
clear all

tic
%acceleration factor:
af = 8;

% our FOV in mm (assuming square FOV):
FOV = 256;

% spiral is inward or outward (inward = 1)
SpiralIn = 0; % only spiral-out is tested and used!

% Do we want to achieve the linear phase ramp with the gradients? (the
% target phase will be close to this, so RF pulse will perform better and
% will need lower max. B1)
LinearPhaseRampWithGrad = 1;

GRT = 10; % gradient raster time in us
dt = 2.5; % RF dwell time in us

% adding some path:
addpath('ESPIRiT');
addpath('ESPIRiT\ESPIRiT_code');
addpath('ESPIRiT\utils');



%% ------------------------------------------------
%  first read in the data from prescan and fixed gradient waveform:
%  ------------------------------------------------

%load('SE-EPI_rawdata.mat');
load('T1wSE_rawdata.mat');
load('gradient_waveform.mat');

% if user wants inward spiral, we do not need the rephasing at the
% beginning (it is at the beginning since we flipped the waveform for inward spiral):
if (SpiralIn==1)
    gx = gx(RampTimes(2)/GRT+1:end);
    gy = gy(RampTimes(2)/GRT+1:end);
end



%% ------------------------------------------------
%  Defining magnetization target and weighting, Calculate target phase
%  ------------------------------------------------

% First, defining the spatial positions, where we want to sample our target

% number of acs lines to use for optimized phase calculation (we want
% it to be a multiple of the acceleration factor)
nacs = round(64/af)*af;

%we want it to be even number:
if (mod(nacs,2)==1)
    nacs = nacs+1;
end

%from this the resolution:
res = FOV/nacs;

% spatial positions vector in read direction (in meter):
x_read = (-FOV/2 : res : FOV/2 - res)'/1000;

% spatial positions vector in phase direction (in meter):
x_phase = (-FOV/2 : res : FOV/2 - res)'/1000;

%long spatial position vector containing all the 2D vectors
% (first coordinate is read direction)
x = zeros(length(x_read)*length(x_phase),2);
x(:,1) = repmat(x_read,[length(x_phase), 1]);

for l=1:length(x_phase)
    x( (l-1)*length(x_read) + 1 : l*length(x_read) , 2)=ones(length(x_read),1) * x_phase(l);
end

% This is the heart of the script: calculating target (including phase-step decreasing)
[phaseTarget, phaseTargetNoSmoothSearch, phaseOpt, myWeights] = optimal_phase_calculation(kspace,af);

%smoothing a bit:
gaussfilter = fspecial('gaussian', 10, 1.5);
target = exp(1i*angle(roifilt2(gaussfilter,real(phaseTarget),myWeights) + 1i*roifilt2(gaussfilter,imag(phaseTarget),myWeights)));

%reshaping target to vector:
target2D = reshape(target,[length(x),1]);

%reshaping signal mask to vector. Now the weighting of pulse design will be
%simple: 1 where there is signal, 0 where there is no signal.
weights_vector = reshape(myWeights,[length(x),1]);

% calculating weight matrix:
W = sparse(diag(weights_vector));

%saving, just in case, and for future use:
save('target.mat','target','myWeights');


%% ------------------------------------------------
%  Defining excitation k-space trajectory
%  ------------------------------------------------

gammabar = 42.577;

% setting ramp times depending on spiral direction:
if (SpiralIn==1)
    RampUpTime = RampTimes(1); % before RF can be turned on we have to ramp up the gradient for inward spiral
    RampDownTime = 0;           % no need for rampdown or rephase in an inward spiral
else
    RampUpTime = 0;             % no need for rampup for outward spiral
    RampDownTime = RampTimes(1)+RampTimes(2);    %ramp down + rephase for k=0;
end

% If user want to create the linear phase ramp with gradients, we delete
% the rephasing gradient blip to k=0, and apply a blip which will take us
% to the appropriate k value for the ramp for spiral out. for spiral in, no
% need to remove rephasing as there is no rephase this case.
if (LinearPhaseRampWithGrad == 1)
    

    % removing rephasing for spiral-out trajectory
    if (SpiralIn == 0)
        gx = gx(1:end-RampTimes(2)/GRT);
        gy = gy(1:end-RampTimes(2)/GRT);
    end

    % to calculate the needed gradient ramp we have to calculate, where are
    % we at the end of the pulse, not taking into account any blip or rephasing.
    % in case of spiral-in, this is zero.

    if (SpiralIn == 0)
        k0(1) = - gammabar * GRT * 1e-3 * sum(gx);
        k0(2) = - gammabar * GRT * 1e-3 * sum(gy);
    else
        k0 = [0, 0];
    end

    %calculate what k-value we need for the linear phaseramp
    kr(1) = - af/4 * 1/FOV * 1e3;
    kr(2) = 0;

    % calculate the needed change in k-space position between end of the pulse
    % and the point according to the linear phase ramp:
    Dk = kr-k0;

    % calculate the gradient blip to achieve this. We want triangular blip with
    % uniform slew rate.
    MaxSlewRate = 188.0; %mT/m/ms

    % how many samples do we need for one side of the triangular blip? (total
    % blip samples will be 2*NumGradBlipSamples + 1 plus a zero at the end ( we
    % do not need a zero at the start because the gradient waveform without the
    % rephasing ends with a zero)
    NumGradBlipSamples = ceil( sqrt(max(abs(Dk))/(gammabar*GRT*1e-6*GRT*MaxSlewRate)) - 1 );

    % what is the exact slew rate with this number of elements for each direction?
    SlewRateUsed = - Dk / (gammabar*GRT*1e-3*(NumGradBlipSamples+1)^2) * (1e3/GRT); % mT/m/ms
    
    % give an error if max. blip amplitude in any direction exceeds max. gradient strength
    % (this would mean we cannot achieve the needed traverse in k-space
    % with a triangular blip. Right now this script cannot handle this
    % case.
    if ( max(abs(( NumGradBlipSamples + 1 ) * SlewRateUsed * GRT * 1e-3 )) > 80.0 )
        disp('Error: Linear phaseramp cannot be achieved with triangular blip!')
        return
    end

    % calculating blip waveform (mT/m/ms):
    for index=1:NumGradBlipSamples+1
        GradBlip(index,1:2) = index*SlewRateUsed*GRT*1e-3;
    end
    for index=NumGradBlipSamples+2:2*NumGradBlipSamples+1
        GradBlip(index,1:2) = GradBlip(NumGradBlipSamples+1,1:2) - ( index - (NumGradBlipSamples+1) )*SlewRateUsed*GRT*1e-3;
    end

    % putting a zero at the end:
    GradBlip = [GradBlip; 0 0];

    % addign this blip to the end of the gradient waveforms:
    gx = [gx; GradBlip(:,1)];
    gy = [gy; GradBlip(:,2)];
    
end


% calculating excitation k-space trajectory, with finer steps than the
% gradient raster time:
for substep = 1:(GRT/dt)
    for index=1:length(gx)
        kx(GRT/dt*(index-1)+substep) = - gammabar * GRT * 1e-3 * ( sum(gx(index+1:end)) + (1-1/(2*GRT/dt)*(2*substep-1))*gx(index) );
        ky(GRT/dt*(index-1)+substep) = - gammabar * GRT * 1e-3 * ( sum(gy(index+1:end)) + (1-1/(2*GRT/dt)*(2*substep-1))*gy(index) );
    end
end

% setting kvalues depending on spiral direction:
if (SpiralIn==1)
    kvalues2D(:,1)=kx(RampUpTime/dt + 1 : end);
    kvalues2D(:,2)=ky(RampUpTime/dt + 1 : end);
else
    kvalues2D(:,1)=kx(1:end-RampDownTime/dt);
    kvalues2D(:,2)=ky(1:end-RampDownTime/dt);
end

% total pulse time in ms:
TotalTime = GRT*length(gx)/1000;


%% -----------------------------------------------------
%  Creating the system matrix of the RF pulse (Fourier matrix)
%  -----------------------------------------------------

A = 1i * dt * exp( -1i * 2*pi * x * (kvalues2D') );

if (SpiralIn==1)
    time = 0:dt:TotalTime*1000-RampUpTime-dt;
else
    time = 0:dt:TotalTime*1000-RampDownTime-dt;
end

%% -----------------------------------------------------
%  set regularization parameters
%  -----------------------------------------------------
beta = dt/10;


lambda = zeros(length(kvalues2D),1);

Matrix1 = A' * W; 
vektor = Matrix1 * target2D;
Matrix2 = Matrix1 * A + eye(length(kvalues2D))*beta;

% loop until abs(sum(b)/max(abs(b))) is large enough, i. e. the needed peak
% B1 amplitude for certain flip angle is low enough to be playable. In each
% step, find and penalize the time points of high B1 amplitudes.

AmplIntOfNormed = 0;
iter=1;

while (iter<31); %this is purely experimental; 31 iteration is usually enough for playable B1 peak.


L = diag(lambda);

lambdaVec(:,iter)=lambda;

%% -----------------------------------------------------
%  Calculating the optimal pulse (using direct inversion)
%  -----------------------------------------------------

Matrix = Matrix2 + L;

% this line is the actual pulse design:
b = linsolve(Matrix,vektor);

AmplIntOfNormed = abs(sum(b)/max(abs(b)));

%storing the pulse integral:
AmplIntOfNormedVec(iter)=AmplIntOfNormed;

% find where the pulse peaks were high:
peak_indeces = find(abs(b)/max(abs(b))>0.9);

%penalize these points:
lambda(peak_indeces) = lambda(peak_indeces) + 0.1;

iter = iter+1;

end

% display the pulse and its results:
figure
subplot(121)
plot(abs(b))
title('pulse magn')

subplot(122)
plot(angle(b))
title('pulse phase')

result = reshape(A*b,[FOV/res,FOV/res]);

figure
subplot(231)
imshow(abs(result),[0 1.2])
title('resulting magnitude map')
colorbar

subplot(232)
imshow((angle(result)),[-pi pi])
title('resulting phase map [rad]')
colormap hsv
colorbar

subplot(234)
imshow(abs(reshape(target2D,[FOV/res,FOV/res])),[0 1.2])
title('target magnitude map')
colorbar

subplot(235)
imshow(angle(reshape(target2D,[FOV/res,FOV/res])).*myWeights,[-pi pi])
title('target phase map []')
colormap hsv
colorbar

subplot(236)
imshow(angle(reshape(target2D,[FOV/res,FOV/res]).*conj(result)).*myWeights,[-0.2 0.2])
title('phase difference, target - result [rad]')
colormap hsv
colorbar
shg;

subplot(233)
imshow(abs(reshape(target2D,[FOV/res,FOV/res]))-abs(result).*myWeights,[-0.1 0.1])
title('magnitude difference, target - result [a.u.]')
colormap hsv
colorbar
shg;

% also plot gradient waveform, for visual check:
figure
subplot(221)
plot(GRT*(1:length(gx)),gx)
title('X gradient pulseshape')
xlabel('time [us]');
ylabel('G_x [mT/m]');
subplot(222)
plot(GRT*(1:length(gx)),gy)
title('Y gradient pulseshape')
xlabel('time [us]');
ylabel('G_y [mT/m]');
subplot(223)
plot(kx,ky)
title('excitation k-space trajectory')
xlabel('k_x [1/m]');
ylabel('k_y [1/m]');
axis equal

% we are done, check calculation time.
toc