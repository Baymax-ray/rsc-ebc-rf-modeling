function [score, details] = scoreRfCv(shapeParams, Egorelative_closest, fired, ...
    trainIdx, testIdx, receptive_field_type, calibration)
%SCORERFCV Profile alpha/offset on train frames and score held-out frames.

fired = fired(:) ~= 0;
nFrames = numel(fired);
trainIdx = toLogicalIndex(trainIdx, nFrames);

if nargin < 5 || isempty(testIdx)
    testIdx = false(nFrames, 1);
else
    testIdx = toLogicalIndex(testIdx, nFrames);
end

rfValue = computeRfValueForClosest(receptive_field_type, Egorelative_closest, shapeParams);
fit = fitLogisticCalibration(rfValue(trainIdx), fired(trainIdx), calibration);

trainEta = fit.alpha * rfValue(trainIdx) - fit.offset;
trainLogLikelihood = logisticLogLikelihood(fired(trainIdx), trainEta);
score = -trainLogLikelihood;

testLogLikelihood = [];
if any(testIdx)
    testEta = fit.alpha * rfValue(testIdx) - fit.offset;
    testLogLikelihood = logisticLogLikelihood(fired(testIdx), testEta);
end

details = struct();
details.rfValue = rfValue;
details.calibrationFit = fit;
details.trainLogLikelihood = trainLogLikelihood;
details.testLogLikelihood = testLogLikelihood;
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
