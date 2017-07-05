% calculating the sum-of-squared VCC g-factors in overlapping voxels with certain background
% phase, used to penalize the high g-factors for minimum search. As used in
% the manuscript.

function  res =  g_factor_cov_fun(phi,C);

nc = size(C,1);

D = C*diag([1 exp(1i*phi)]);

D(1+nc:2*nc,:) = conj(D(1:nc,:));


% penalize the sum of g-factor squares

res = trace(pinv(D'*D).*(D'*D));


