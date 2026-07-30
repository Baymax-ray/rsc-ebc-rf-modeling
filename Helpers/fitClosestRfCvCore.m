function result = fitClosestRfCvCore(Egorelative_closest, fired, receptive_field_type, ...
    numClosest, calibration, folds, foldNaiveLogLikelihood, fitOptions, inheritResult)
%FITCLOSESTRFCVCORE Fit one calibration, optionally inheriting valid folds.

if nargin < 8 || isempty(fitOptions)
    fitOptions = struct();
end
if nargin < 9
    inheritResult = [];
end

fired = fired(:) ~= 0;
nFrames = numel(fired);
spec = getRfSpec(receptive_field_type);
foldNaiveLogLikelihood = foldNaiveLogLikelihood(:);

if numel(folds) ~= numel(foldNaiveLogLikelihood)
    error('folds and foldNaiveLogLikelihood must have matching lengths.');
end

foldParams = zeros(numel(folds), spec.nParams);
foldAlpha = zeros(numel(folds), 1);
foldOffset = zeros(numel(folds), 1);
foldLogLikelihood = zeros(numel(folds), 1);
foldTrainLogLikelihood = zeros(numel(folds), 1);
foldCalibrationMethod = strings(numel(folds), 1);

for k = 1:numel(folds)
    if canInheritFold(calibration, inheritResult, k)
        fprintf('  Fold %d/%d (%s inherited from glm_free)\n', k, numel(folds), calibration);
        foldParams(k, :) = inheritResult.foldParams(k, :);
        foldAlpha(k) = inheritResult.foldAlpha(k);
        foldOffset(k) = inheritResult.foldOffset(k);
        foldLogLikelihood(k) = inheritResult.foldLogLikelihood(k);
        foldTrainLogLikelihood(k) = inheritResult.foldTrainLogLikelihood(k);
        foldCalibrationMethod(k) = "inherited_from_glm_free";
        continue
    end

    fprintf('  Fold %d/%d (%s)\n', k, numel(folds), calibration);
    setDeterministicSeed(fitOptions, receptive_field_type, sprintf('fold_%d', k));
    [shapeParams, ~] = fitRfShapeOnTrain(Egorelative_closest, fired, ...
        folds(k).trainIdx, receptive_field_type, calibration, fitOptions);
    [~, details] = scoreRfCv(shapeParams, Egorelative_closest, fired, ...
        folds(k).trainIdx, folds(k).testIdx, receptive_field_type, calibration);

    foldParams(k, :) = shapeParams;
    foldAlpha(k) = details.calibrationFit.alpha;
    foldOffset(k) = details.calibrationFit.offset;
    foldLogLikelihood(k) = details.testLogLikelihood;
    foldTrainLogLikelihood(k) = details.trainLogLikelihood;
    foldCalibrationMethod(k) = string(details.calibrationFit.method);
end

cvLogLikelihood = sum(foldLogLikelihood);
naiveCvLogLikelihood = sum(foldNaiveLogLikelihood);

fullFitParams = [];
fullFitAlpha = [];
fullFitOffset = [];
fullFitLogLikelihood = NaN;
AIC = NaN;
BIC = NaN;
params = [];

if getOption(fitOptions, 'DoFullFit', true)
    if canInheritFull(calibration, inheritResult)
        fprintf('  Full fit (%s inherited from glm_free)\n', calibration);
        fullFitParams = inheritResult.fullFitParams;
        fullFitAlpha = inheritResult.fullFitAlpha;
        fullFitOffset = inheritResult.fullFitOffset;
        fullFitLogLikelihood = inheritResult.loglikelihood;
        AIC = inheritResult.AIC;
        BIC = inheritResult.BIC;
        params = inheritResult.params;
    else
        allTrainIdx = true(nFrames, 1);
        fprintf('  Full fit (%s)\n', calibration);
        setDeterministicSeed(fitOptions, receptive_field_type, 'full');
        [fullFitParams, ~] = fitRfShapeOnTrain(Egorelative_closest, fired, ...
            allTrainIdx, receptive_field_type, calibration, fitOptions);
        [~, fullDetails] = scoreRfCv(fullFitParams, Egorelative_closest, fired, ...
            allTrainIdx, [], receptive_field_type, calibration);
        fullFitAlpha = fullDetails.calibrationFit.alpha;
        fullFitOffset = fullDetails.calibrationFit.offset;
        fullFitLogLikelihood = fullDetails.trainLogLikelihood;
        params = [fullFitParams, fullFitAlpha, fullFitOffset];

        nModelParams = spec.nParams + 2;
        AIC = 2 * nModelParams - 2 * fullFitLogLikelihood;
        BIC = log(nFrames) * nModelParams - 2 * fullFitLogLikelihood;
    end
end

result = emptyCvResultStruct(1);
result.frame = nFrames;
result.loglikelihood = fullFitLogLikelihood;
result.receptive_field_type = receptive_field_type;
result.AIC = AIC;
result.BIC = BIC;
result.params = params;
result.reFireRate = [];
result.pvalue = [];
result.cvLogLikelihood = cvLogLikelihood;
result.cvLogLikelihoodPerFrame = cvLogLikelihood / nFrames;
result.naiveCvLogLikelihood = naiveCvLogLikelihood;
result.cvDeltaLogLikelihood = cvLogLikelihood - naiveCvLogLikelihood;
result.foldLogLikelihood = foldLogLikelihood;
result.foldNaiveLogLikelihood = foldNaiveLogLikelihood;
result.foldTrainLogLikelihood = foldTrainLogLikelihood;
result.foldParams = foldParams;
result.foldAlpha = foldAlpha;
result.foldOffset = foldOffset;
result.foldCalibrationMethod = foldCalibrationMethod;
result.calibration = calibration;
result.numClosest = numClosest;
result.fullFitParams = fullFitParams;
result.fullFitAlpha = fullFitAlpha;
result.fullFitOffset = fullFitOffset;
end

function tf = canInheritFold(calibration, inheritResult, k)
tf = strcmp(char(string(calibration)), 'glm_pos') && ~isempty(inheritResult) && ...
    numel(inheritResult.foldAlpha) >= k && isfinite(inheritResult.foldAlpha(k)) && ...
    inheritResult.foldAlpha(k) >= 0;
end

function tf = canInheritFull(calibration, inheritResult)
tf = strcmp(char(string(calibration)), 'glm_pos') && ~isempty(inheritResult) && ...
    ~isempty(inheritResult.fullFitAlpha) && isfinite(inheritResult.fullFitAlpha) && ...
    inheritResult.fullFitAlpha >= 0;
end

function setDeterministicSeed(fitOptions, receptive_field_type, foldIdentity)
if isstruct(fitOptions) && isfield(fitOptions, 'SeedTaskId') && ...
        isfield(fitOptions, 'SeedNeuronName')
    seed = makeDeterministicSeed(fitOptions.SeedTaskId, ...
        fitOptions.SeedNeuronName, receptive_field_type, foldIdentity);
    rng(seed, "twister");
end
end

function value = getOption(options, name, defaultValue)
if isstruct(options) && isfield(options, name) && ~isempty(options.(name))
    value = options.(name);
else
    value = defaultValue;
end
end
