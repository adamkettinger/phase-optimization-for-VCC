%this is a sample script to compare GRAPPA and VCC-GRAPPA reconstructions and g-factors of
% conventional and optimized scans.

% note that this is a very simple example with non-iterative GRAPPA, only
% to show some insight to VCC reconstruction.
% therefore there are severe remaining aliasing artefacts, and g-factors
% are different from what they would be with an iterative and artefact-free
% recon (as they are in the manuscript). The iterative GRAPPA code 
% used in the manuscript could not be made public due to copyright.
% However, the main effect (some benefit using VCC alone,
% more benefit using VCC on the tailored scan) is clearly ovservable,
% and the reconstruction could be extended to iterative GRAPPA.

% Written by Adam Kettinger, 2017. Budapest, Hungary

clear all

% load data:
load('T1wSE_rawdata.mat');
%load('SE-EPI_rawdata.mat');

% set acceleration factor and ACS line number:
af = 8;
nacs = 32;

[S1 S2 S3] = size(kspace_conv);

% simulate accelerated k-space:
kspace_conv_folded = kspace_conv(:,1:af:end,:);
kspace_opt_folded = kspace_opt(:,1:af:end,:);

% simulate acs datasets:
kspace_conv_acs = kspace_conv(:,S2/2+1-nacs/2:S2/2+nacs/2,:);
kspace_opt_acs = kspace_opt(:,S2/2+1-nacs/2:S2/2+nacs/2,:);

% create VCC signals:
kspace_vcc_conv = VCC_signal_creation( kspace_conv );
kspace_vcc_opt  = VCC_signal_creation( kspace_opt );

%simulating folded and acs datasets for VCC signals:
kspace_vcc_conv_folded = VCC_signal_creation(kspace_conv_folded);
kspace_vcc_opt_folded = VCC_signal_creation(kspace_opt_folded);
kspace_vcc_conv_acs = VCC_signal_creation(kspace_conv_acs);
kspace_vcc_opt_acs = VCC_signal_creation(kspace_opt_acs);


% do the reconstructions. first, standard and VCC recon on conventional
% measurement:
[recon_conv , gfact_conv ] = grappa(kspace_conv_folded,kspace_conv_acs,af);
[recon_conv_vcc , gfact_conv_vcc ] = grappa(kspace_vcc_conv_folded,kspace_vcc_conv_acs,af);

%now, standard and VCC recon on tailored measurement:
[recon_opt , gfact_opt ] = grappa(kspace_opt_folded,kspace_opt_acs,af);
[recon_opt_vcc , gfact_opt_vcc ] = grappa(kspace_vcc_opt_folded,kspace_vcc_opt_acs,af);

% displaying results. take into accound that VCC reconstructions have
% higher intensity due to 2x number of channels:

figure
subplot(221)
imshow(recon_conv, [0 6e-6]);
title('standard recon, conventional scan');
subplot(222)
imshow(recon_conv_vcc/sqrt(2), [0 6e-6]);
title('VCC recon, conventional scan');
subplot(223)
imshow(recon_opt, [0 6e-6]);
title('standard recon, tailored scan');
subplot(224)
imshow(recon_opt_vcc/sqrt(2), [0 6e-6]);
title('VCC recon, tailored scan');

figure
subplot(221)
imshow(gfact_conv, [1 7]);
title('gfactors of standard recon, conventional scan');
colorbar
subplot(222)
imshow(gfact_conv_vcc/sqrt(2), [1 7]);
title('gfactors of VCC recon, conventional scan');
colorbar
subplot(223)
imshow(gfact_opt, [1 7]);
title('gfactors of standard recon, tailored scan');
colorbar
subplot(224)
imshow(gfact_opt_vcc/sqrt(2), [1 7]);
title('gfactors of VCC recon, tailored scan');
colorbar
colormap jet

figure
subplot(221)
imshow(gfact_conv, [1 3]);
title('gfactors of standard recon, conventional scan');
colorbar
subplot(222)
imshow(gfact_conv_vcc/sqrt(2), [1 3]);
title('gfactors of VCC recon, conventional scan');
colorbar
subplot(223)
imshow(gfact_opt, [1 3]);
title('gfactors of standard recon, tailored scan');
colorbar
subplot(224)
imshow(gfact_opt_vcc/sqrt(2), [1 3]);
title('gfactors of VCC recon, tailored scan');
colorbar
colormap jet

