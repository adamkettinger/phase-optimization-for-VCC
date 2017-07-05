function [combImage,cMap,ESP,weights] = cMapEspirit(data,acsSize);

% Computes Coil Sensitivity Maps using ESPIRiT with code from Martin
% Uecker and Micheal Lustig
%
% input:    'data'      k-space data (dimensions: Nread,Nphase,Ncoil);
%           'acsSize'   size of calibration region (e.g. [64, 32]);
%
% output:   'combImage' combined image
%           'cMap'      coil sensitivities
%           'ESP'       ESPIRiT operation
%
% 27 Nov 2015, M. Blaimer
%

    % ESPIRIT variables
    ksize = [5,5];                  % Kernel size

    % Threshold for picking singular vercors of the calibration matrix
    % (relative to largest singlular value).
    eigThresh_1 = 0.02;

    % Threshold of eigen vector decomposition in image space. 
    eigThresh_im = 0.98;

    calib = crop(data,[acsSize(1),acsSize(2),size(data,3)]);   % extract calibration data
    [k,S] = dat2Kernel(calib,ksize);                           % calibration
    idx = max(find(S >= S(1)*eigThresh_1));                    % check eigenvalues
    [M,W] = kernelEig(k(:,:,:,1:idx),[size(data,1), size(data,2)]);
 
    %% Compute ESPIRiT Maps 

    % Crop sensitivity maps according to eigenvalues==1
    cMap = M(:,:,:,end);


    % Weight the eigenvectors with soft-senses eigen-values
    weights = W(:,:,end) ;
    weights = (weights - eigThresh_im)./(1-eigThresh_im).* (W(:,:,end) > eigThresh_im);
    weights = -cos(pi*weights)/2 + 1/2;


    % Create ESPIRiT operator
    ESP = ESPIRiT(cMap,weights);

    %% Coil combination
    dataIm = data;
    dataIm = ifftshift(ifft(ifftshift(dataIm,1),[],1),1);       % iFFT along read
    dataIm = ifftshift(ifft(ifftshift(dataIm,2),[],2),2);       % iFFT along PE

    % combined Image
    combImage = ESP'*dataIm;
    
    
   
end % function