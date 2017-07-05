% penalty function for searching smooth optimal target, as used in
% manuscript.

function res = PhaseStepsMinimize(phaseShift, phaseOpt, myWeights, af)


% assuming both phaseTarget and phaseOpt are real (angles).

% extend phase shift to all PE FOV:
phaseShiftExtended = repmat(phaseShift,[af,1]);


% calculate phase gradients (penalty function is square of phase gradient
% magnitude). take into account that i.e. 2*pi - epsilon gradient is
% epsilon gradient!
phasegrad1 = myWeights(2:end,:).*diff((phaseShiftExtended + phaseOpt), 1, 1);
phasegrad1 = angle(exp(1i*phasegrad1));

% same with second direction
phasegrad2 = myWeights(:,2:end).*diff((phaseShiftExtended + phaseOpt), 1, 2);
phasegrad2 = angle(exp(1i*phasegrad2));

%penalize high phase gradients:
res = sum(sum(phasegrad1.^2)) + sum(sum(phasegrad2.^2));
