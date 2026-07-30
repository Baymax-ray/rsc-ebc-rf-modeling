function [foldNaiveLogLikelihood, naiveCvLogLikelihood] = computeFoldNaiveLogLikelihood(fired, folds)
%COMPUTEFOLDNAIVELOGLIKELIHOOD Score train-rate baselines for fixed folds.

fired = fired(:) ~= 0;
nFrames = numel(fired);
foldNaiveLogLikelihood = zeros(numel(folds), 1);

for k = 1:numel(folds)
    trainIdx = toLogicalIndex(folds(k).trainIdx, nFrames);
    testIdx = toLogicalIndex(folds(k).testIdx, nFrames);
    trainRate = mean(fired(trainIdx));
    naiveEta = logitClipped(trainRate) * ones(sum(testIdx), 1);
    foldNaiveLogLikelihood(k) = logisticLogLikelihood(fired(testIdx), naiveEta);
end

naiveCvLogLikelihood = sum(foldNaiveLogLikelihood);
end

function idx = toLogicalIndex(idx, nFrames)
if islogical(idx)
    idx = idx(:);
else
    tmp = false(nFrames, 1);
    tmp(idx(:)) = true;
    idx = tmp;
end

if numel(idx) ~= nFrames
    error('Index length does not match number of frames.');
end
end

function x = logitClipped(p)
epsval = 1e-6;
p = min(max(p, epsval), 1 - epsval);
x = log(p / (1 - p));
end
