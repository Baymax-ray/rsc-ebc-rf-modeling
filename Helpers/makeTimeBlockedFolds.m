function folds = makeTimeBlockedFolds(nFrames, nFolds)
%MAKETIMEBLOCKEDFOLDS Build contiguous time-block cross-validation folds.
%
% folds(k).testIdx is one contiguous block. folds(k).trainIdx is every
% frame outside that block.

if nargin < 2 || isempty(nFolds)
    nFolds = 5;
end

nFrames = double(nFrames);
nFolds = double(nFolds);

if ~isscalar(nFrames) || nFrames < 2 || fix(nFrames) ~= nFrames
    error('nFrames must be an integer >= 2.');
end
if ~isscalar(nFolds) || nFolds < 2 || fix(nFolds) ~= nFolds
    error('nFolds must be an integer >= 2.');
end

nFolds = min(nFolds, nFrames);
baseSize = floor(nFrames / nFolds);
foldSizes = baseSize * ones(1, nFolds);
foldSizes(1:mod(nFrames, nFolds)) = foldSizes(1:mod(nFrames, nFolds)) + 1;

folds = repmat(struct('fold', [], 'trainIdx', [], 'testIdx', [], ...
    'testStart', [], 'testEnd', []), nFolds, 1);

testStart = 1;
for k = 1:nFolds
    testEnd = testStart + foldSizes(k) - 1;
    testIdx = false(nFrames, 1);
    testIdx(testStart:testEnd) = true;

    folds(k).fold = k;
    folds(k).testIdx = testIdx;
    folds(k).trainIdx = ~testIdx;
    folds(k).testStart = testStart;
    folds(k).testEnd = testEnd;

    testStart = testEnd + 1;
end
end
