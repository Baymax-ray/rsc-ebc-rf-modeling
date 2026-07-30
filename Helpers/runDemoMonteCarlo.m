function runDemoMonteCarlo(demoData, config)
%RUNDEMOMONTECARLO Compute CV null p-values and regenerated firing rates.

for rfIdx = 1:numel(config.rfTypes)
    rfType = config.rfTypes{rfIdx};
    for calIdx = 1:numel(config.calibrations)
        calibration = config.calibrations{calIdx};
        groupName = sprintf('%s_%dclosest_cv_%s', ...
            rfType, config.numClosest, calibration);
        bestFile = fullfile(config.resultDir, [groupName '_best.mat']);
        if ~isfile(bestFile)
            error('Best result not found: %s', bestFile);
        end

        S = load(bestFile, 'results_best');
        results_best = S.results_best;
        if numel(results_best) ~= 1
            error('The Demo expects exactly one result in %s.', bestFile);
        end
        entry = results_best(1);
        validateEntry(entry, demoData, rfType, calibration);

        [oofProbability, reconstructedFoldLL] = ...
            reconstructOofProbability(entry, demoData, rfType);
        validateReconstructedCv(entry, reconstructedFoldLL);
        nullLogLikelihood = computeNullLogLikelihoods( ...
            oofProbability, demoData.fired, demoData.folds, ...
            config.nPermutations, config.useParallel, ...
            config.monteCarloSeedBase);
        pValue = sum(nullLogLikelihood >= entry.cvLogLikelihood) / ...
            config.nPermutations;

        rfValue = computeRfValueForClosest( ...
            rfType, demoData.egoClosest, entry.fullFitParams);
        eta = entry.fullFitAlpha * rfValue - entry.fullFitOffset;
        firingProbability = clippedLogistic(eta);
        stream = RandStream('mt19937ar', ...
            'Seed', config.regeneratedSeed);
        regenerated = rand(stream, numel(firingProbability), 1) < ...
            firingProbability;

        results_best(1).pvalue = pValue;
        results_best(1).reFireRate = mean(regenerated) * config.fps;
        filledFile = fullfile(config.resultDir, [groupName '_filled.mat']);
        save(filledFile, 'results_best');

        fprintf('%s: CV LL %.3f, p %.4g, regenerated %.3f Hz\n', ...
            groupName, entry.cvLogLikelihood, pValue, ...
            results_best(1).reFireRate);
    end
end
end

function validateEntry(entry, demoData, rfType, calibration)
required = {'name', 'frame', 'cvLogLikelihood', 'foldLogLikelihood', ...
    'foldParams', 'foldAlpha', 'foldOffset', 'fullFitParams', ...
    'fullFitAlpha', 'fullFitOffset', 'calibration', ...
    'receptive_field_type'};
for i = 1:numel(required)
    if ~isfield(entry, required{i}) || isempty(entry.(required{i}))
        error('Best result is missing required field %s.', required{i});
    end
end
if ~strcmp(char(string(entry.name)), demoData.name)
    error('Result neuron name does not match the Demo data.');
end
if entry.frame ~= demoData.frame
    error('Result frame count does not match the Demo data.');
end
if ~strcmp(char(string(entry.calibration)), calibration)
    error('Calibration mismatch: expected %s.', calibration);
end
if ~strcmp(char(string(entry.receptive_field_type)), rfType)
    error('RF type mismatch: expected %s.', rfType);
end
end

function [oofProbability, foldLogLikelihood] = ...
    reconstructOofProbability(entry, demoData, rfType)
nFrames = demoData.frame;
nFolds = numel(demoData.folds);
if numel(entry.foldLogLikelihood) ~= nFolds || ...
        size(entry.foldParams, 1) ~= nFolds || ...
        numel(entry.foldAlpha) ~= nFolds || ...
        numel(entry.foldOffset) ~= nFolds
    error('Stored CV fields do not contain %d folds.', nFolds);
end

fired = demoData.fired(:) ~= 0;
oofProbability = nan(nFrames, 1);
foldLogLikelihood = zeros(nFolds, 1);
for k = 1:nFolds
    rfValue = computeRfValueForClosest( ...
        rfType, demoData.egoClosest, entry.foldParams(k, :));
    testIdx = demoData.folds(k).testIdx;
    eta = entry.foldAlpha(k) * rfValue(testIdx) - entry.foldOffset(k);
    probability = clippedLogistic(eta);
    oofProbability(testIdx) = probability;
    foldLogLikelihood(k) = logisticLogLikelihood( ...
        fired(testIdx), eta);
end
if any(~isfinite(oofProbability))
    error('Failed to reconstruct finite OOF probabilities for every frame.');
end
end

function validateReconstructedCv(entry, reconstructedFoldLL)
storedFoldLL = entry.foldLogLikelihood(:);
foldTolerance = 1e-8 * max(1, abs(storedFoldLL));
if any(abs(reconstructedFoldLL - storedFoldLL) > foldTolerance)
    error('Reconstructed fold LL does not match foldLogLikelihood.');
end
totalTolerance = 1e-8 * max(1, abs(entry.cvLogLikelihood));
if abs(sum(storedFoldLL) - entry.cvLogLikelihood) > totalTolerance || ...
        abs(sum(reconstructedFoldLL) - entry.cvLogLikelihood) > totalTolerance
    error('Reconstructed CV LL does not match cvLogLikelihood.');
end
end

function nullLogLikelihood = computeNullLogLikelihoods( ...
    oofProbability, fired, folds, nPermutations, useParallel, seedBase)
nFolds = numel(folds);
testIndices = cell(nFolds, 1);
trainRates = zeros(nFolds, 1);
fired = fired(:) ~= 0;
for k = 1:nFolds
    testIndices{k} = find(folds(k).testIdx);
    trainRates(k) = mean(fired(folds(k).trainIdx));
end

nullLogLikelihood = zeros(nPermutations, 1);
if useParallel
    parfor i = 1:nPermutations
        nullLogLikelihood(i) = oneNullLogLikelihood( ...
            oofProbability, testIndices, trainRates, seedBase + i);
    end
else
    for i = 1:nPermutations
        nullLogLikelihood(i) = oneNullLogLikelihood( ...
            oofProbability, testIndices, trainRates, seedBase + i);
    end
end
end

function value = oneNullLogLikelihood( ...
    oofProbability, testIndices, trainRates, seed)
stream = RandStream('mt19937ar', 'Seed', seed);
value = 0;
for k = 1:numel(testIndices)
    idx = testIndices{k};
    randomFired = rand(stream, numel(idx), 1) < trainRates(k);
    probability = oofProbability(idx);
    value = value + sum(randomFired .* log(probability) + ...
        (1 - randomFired) .* log(1 - probability));
end
end

function probability = clippedLogistic(eta)
probability = 1 ./ (1 + exp(-eta));
epsValue = 1e-15;
probability = min(max(probability, epsValue), 1 - epsValue);
end
