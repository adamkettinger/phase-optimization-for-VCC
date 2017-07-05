function [ kspace_vcc ] = VCC_signal_creation( kspace )
% script for creating virtual coil signals, that could be used with any
% further reconstruction, e.g. VCC-GRAPPA.
%
% input: a matrix with k-space data of one slice, with dimension order:
% coil, PE, RO
%
% output: a matrix containing the original k-space data as well as virtual
% coil data. virtual coil signals are concatenated to the first dimenision,
% i. e. the output will have 2xNc coils. Virtual coil signals are stored
% with the same order as the originals, i. e. coil with index Nc+1 will be
% the corresponding virtual coil of the physical coil with index 1.

% written by Adam Kettinger, 2017. Budapest, Hungary

[nC, nPE, nRO] = size(kspace);

VCC_signals = conj(flip(flip(kspace,2),3));

% take care of k-space center!
if mod(nPE,2) == 0
    VCC_signals = circshift(VCC_signals,[0 1 0]);
end
if mod(nRO,2) == 0
    VCC_signals = circshift(VCC_signals,[0 0 1]);
end

kspace_vcc = [kspace; VCC_signals];

end

