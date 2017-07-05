function [phaseTarget2, phaseTarget, phaseOpt, myWeights] = optimal_phase_calculation(kspace,af)
% This script computes a target phase pattern
% for opimized phase-constrained parallel MRI reconstructions, as used in
% the manuscript.
% Ideally, when applying an acceleration factor R the g-factor of the optimized
% reconstruction will be as low as for an R/2 acceleration
%
% original script M.Blaimer, 09. Mar 2016
%
% input:    - matrix of k-space
%           - acceleration factor
% 
% output:   - phaseTarget2: target phase map
%           - phaseTarget: target phase map without searhc for smoother equivalent (for testing)
%           - phaseOpt: Optimal phase map (target + intrinsic object phase)
%           - myWeights: mask of original image (1 for relevant signal)


% ----------
% Initialize
%

display('Preparation ...');

% calculate number of lines in phase optimization:
nacs = round(64/af)*af;
%we want it to be even number:
if (mod(nacs,2)==1)
    nacs = nacs+1;
end

% --------------------------------------------------------------------------
% Preparation

kspace      = ifftshift(ifft(ifftshift(kspace,2),[],2),2);  % remove read oversampling
S2          = size(kspace,2);
kspace      = kspace(:,S2/4+1:S2*3/4,:);                  % remove read oversampling
kspace      = fftshift(fft(fftshift(kspace,2),[],2),2);     % remove read oversampling

data        = permute(kspace,[3 1 2]);                      % bring data in order: channels, lines, readout, slices

% extract fully sampled ACS region with low resolution (e.g. 32 x 32)
% this data will be used for computing the phase pattern)
ny_temp = size(data,2);
nx_temp = size(data,3);
acs = data(:,ny_temp/2+1-nacs/2:ny_temp/2+nacs/2,nx_temp/2+1-nacs/2:nx_temp/2+nacs/2);

nc = size(acs,1)                  % number of coils
ny = size(acs,2)                  % number of PE lines
nx = size(acs,3)                  % number of readout points

% ------------------------------ Start: prepare images --------------------
%
% the goal is to extract coil sensitivity and object phase information

% filter ACS data to avoid Gibbs-Ringing
hannFilter                      = tukeywin(size(acs,2),1)';
xFilter                        = tukeywin(nx,1);

h = repmat(xFilter*hannFilter,[1 1 nc]);
h = permute(h,[3 2 1]);

acsImage = acs.*h;

% Compute coil sensitivity maps using ESPIRiT
[adaptImage,cMap,ESP,cWeights]   = cMapEspirit(permute(acsImage,[3 2 1]),[32 32]);

cMap                        = cMap.*repmat(cWeights,[1 1 nc]);
adaptImage                  = permute(adaptImage,[2 1]);
cMap                        = permute(cMap,[3 2 1]);
cWeights                    = permute(cWeights,[2 1]);

% Extract object phase
phaseImage          = angle(adaptImage);

phaseImage = repmat(phaseImage,[1 1 nc]);
phaseImage = permute(phaseImage,[3 1 2]);
    

% ----------------------- Start: compute optimal phase pattern ------------
%
% the goal is to find an object phase pattern that minimizes the g-Factor
%
% in total, R-1 phase differences delta_phi have to be found

delta_phi   = [1:af-1]*pi/2*1;       % start values (assuming near optimal phase difference of pi/2 between pixels)
phaseOpt    = zeros(ny,nx);            % allocate memory
options=optimset('MaxFunEvals',1e6,'MaxIter',1e6);

%mask where there is signal (we only care about the phase there):
myWeights_raw = logical(abs(adaptImage)>max(max(abs(adaptImage)))*0.07);

%we do not want accidental holes in myWeights:
myWeights = imfill(myWeights_raw,'holes');

%repeat the mask for all coils:
myWeightsRep = permute(repmat(myWeights,1,1,nc),[3 1 2]);

%visual check on the weights:
figure
imshow(myWeights)

%calculating the optimum phase:
for y=1:ny/af,
    
    yy = y:ny/af:ny;    % indices in PE direction
  
        parfor xx = 1:nx
            
            
            if abs(sum(adaptImage(yy,xx),1))>0,

                C = cMap(1:nc,yy,xx).*myWeightsRep(1:nc,yy,xx); % masked coil sensitivity information

                % now we simply minimize sum of g-factor squares in the overlapping voxel group, do not
                % try to make them equal to the R/2 case.
                p = fminsearch(@(p) g_factor_cov_fun(p,C),delta_phi,options);
                
                phi = [0, p];
            
                phaseOpt(yy,xx) = phi+y*pi/2/(ny/af);
                
            else
                % if there is no signal, just use the pi/2 linear phaseramp
                phaseOpt(yy,xx) = y*pi/2/(ny/af)+[0 delta_phi];
            end
        
        end
        

end

% set phaseOpt between -pi and pi:
phaseOpt = angle(exp(1i*phaseOpt));

        %visual check after normalizing phase to -pi..pi
        figure(3), 
        imagesc((phaseOpt)); 
        axis off image; 
        colormap('hot');
        title('optimal pahse distribution'); 
        colorbar;
        drawnow

% target phase pattern for RF pulse: = optimal phase - object phase
se = strel('disk',3);
cWeights = imdilate(cWeights,se);
phaseTarget = cWeights.*(phaseOpt - squeeze(phaseImage(1,:,:)));
phaseTarget = angle(exp(1i*phaseTarget));

% calculating the intra-mask phaseshifts (maintaining the phase difference
% between overlapping voxels, but searching for smoother target):
phaseShift = zeros(size(phaseOpt,1)/af,size(phaseOpt,2));
options = optimset('MaxFunEvals',1e8,'MaxIter',1e8);
phaseShift = fminunc(@(phaseShift) PhaseStepsMinimize(phaseShift,phaseTarget,myWeights,af),zeros(size(phaseTarget,1)/af,size(phaseTarget,2)),options);

% apply the computed phase shifts to target:
phaseTarget2 = cWeights.*angle(exp(1i*(repmat(phaseShift,[af,1])+phaseTarget)));

% % Restrict the mask of sensitivities to the mask of signal:
for k=1:size(cMap,1)
     cMap(k,:,:) = squeeze(cMap(k,:,:)).*myWeights;
end

% Calculate simulated VCC-SENSE g-factors:

% generate effective coil sensitivities with object phasev
aMap                = cMap(1:nc,:,:).*exp(1i.*phaseImage);
aMap(1+nc:2*nc,:,:) = conj(aMap(1:nc,:,:));     % generate virtual conjugate coils
    
% generate effective coil sensitivities with optimized phase
for l=1:nc, bMap(l,:,:) = squeeze(aMap(l,:,:)).*exp(1i.*phaseTarget2); end
bMap(1+nc:2*nc,:,:) = conj(bMap(1:nc,:,:));               % generate virtual conjugate coils

% generate effective coil sensitivities with smoothed optimized phase (RF
% target will be smoothed the same way)
phaseTargetSmoothed = angle(imgaussfilt(cos(phaseTarget2),2)+1i*imgaussfilt(sin(phaseTarget2),2));
for l=1:nc, dMap(l,:,:) = squeeze(aMap(l,:,:)).*exp(1i.*phaseTargetSmoothed); end
dMap(1+nc:2*nc,:,:) = conj(dMap(1:nc,:,:));               % generate virtual conjugate coils

% generate effective coil sensitivities for simple linear phase:
linearPhase = repmat([0:pi*af/2/size(phaseTarget2,1):pi*af/2-pi*af/2/size(phaseTarget2,1)]',[1,size(phaseTarget2,2)]);
for l=1:nc, eMap(l,:,:) = squeeze(aMap(l,:,:)).*exp(1i.*linearPhase); end
eMap(1+nc:2*nc,:,:) = conj(eMap(1:nc,:,:));               % generate virtual conjugate coils


% compute g-Factor maps
g_PC    = gfactor(aMap,af);                   % standard phase-constrained´
g_opt   = gfactor(bMap,af);                   % optimized phase-constrained
g_conv  = gfactor(cMap(1:nc,:,:),af);         % conventional SENSE for comparison
g_conv_half = gfactor(cMap(1:nc,:,:),af/2);   % conventional SENSE with reduced acceleration factor (af/2)
g_smoothed_opt = gfactor(dMap,af);            % optimized phase-constrained with smoothed optimized phase
g_linear = gfactor(eMap,af);                  % phase-constrained with simple linear phase added to object phase


% ---------------------------- Display ------------------------------------
%
figure(1),
imagesc(([abs(g_conv_half) abs(g_conv) abs(g_PC) ; abs(g_opt) abs(g_smoothed_opt) abs(g_linear)]),[0 af]); 
axis image off
colorbar
colormap('jet');
title({['Conv. SENSE (R=',num2str(af/2),')   |   Conv. SENSE (R=',num2str(af),')   |   PC-SENSE (R=',num2str(af),')'] ; ['Opt. PC-SENSE (R=',num2str(af),')   |   Smoothed Opt. PC-SENSE (R=',num2str(af),')   |   Linear phase PC-SENSE (R=',num2str(af),')']});
drawnow;

% g-factors of R/2 and optimized PC (R) are approximately the same?
disp('sum of squared difference between standard R/2 g-factors and optimized PC R g-factors:')
sum(sum((g_opt - g_conv_half).^2))

figure(2), 
imagesc((phaseTarget2)); 
axis off image; 
colormap('hot');
title(y); 
colorbar;
title('Target Phase Pattern');
drawnow

figure
subplot(121)
imshow(phaseTarget.*myWeights,[-pi pi])
title('phase target without phasestep optimization')
subplot(122)
imshow(phaseTarget2.*myWeights,[-pi pi])
title('phase target after phasestep optimization (maintaining g-factors)')
colormap hsv

% make target complex for further use:
phaseTarget2 = exp(1i*phaseTarget2);

figure
imshow(angle(adaptImage),[-pi pi])
colormap hsv
title('phase of combined image [rad]')

figure
imshow(abs(adaptImage),[])
title('magnitude of combined image [a.u.]')


end