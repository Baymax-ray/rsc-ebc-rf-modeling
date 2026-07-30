function selectDemoBestFits(config)
%SELECTDEMOBESTFITS Build fold-stitched best results for configured runs.

for rfIdx = 1:numel(config.rfTypes)
    rfType = config.rfTypes{rfIdx};
    for calIdx = 1:numel(config.calibrations)
        expectedCalibration = config.calibrations{calIdx};
        groupName = sprintf('%s_%dclosest_cv_%s', ...
            rfType, config.numClosest, expectedCalibration);
        groupFiles = configuredGroupFiles( ...
            config.resultDir, groupName, expectedCalibration, config.runIds);

        fprintf('Selecting best fits for %s (%d runs)\n', ...
            groupName, numel(groupFiles));
        data = loadGroupFiles(config.resultDir, groupFiles);
        results_best = buildBestResults( ...
            data, groupFiles, expectedCalibration);

        outFile = fullfile(config.resultDir, [groupName '_best.mat']);
        save(outFile, 'results_best');
    end
end
end

function groupFiles = configuredGroupFiles( ...
    resultDir, groupName, calibration, runIds)
groupFiles = struct('name', {}, 'group', {}, ...
    'calibration', {}, 'suffix', {});
runIds = sort(runIds);

for i = 1:numel(runIds)
    runId = runIds(i);
    fileName = sprintf('%s_%d.mat', groupName, runId);
    if ~isfile(fullfile(resultDir, fileName))
        error('Configured fit result not found: %s', fileName);
    end
    groupFiles(i).name = fileName;
    groupFiles(i).group = groupName;
    groupFiles(i).calibration = calibration;
    groupFiles(i).suffix = runId;
end
end

function data = loadGroupFiles(resultDir, groupFiles)
data = struct('fileName', {}, 'suffix', {}, 'results', {});

for f = 1:numel(groupFiles)
    inFile = fullfile(resultDir, groupFiles(f).name);
    S = load(inFile);
    if ~isfield(S, 'results_to_save')
        error('File %s does not contain results_to_save.', groupFiles(f).name);
    end

    data(f).fileName = groupFiles(f).name;
    data(f).suffix = groupFiles(f).suffix;
    data(f).results = S.results_to_save;
    validateRequiredFields(data(f).results, groupFiles(f).name);
end
end

function results_best = buildBestResults(data, groupFiles, expectedCalibration)
baseResults = data(1).results;
results_best = baseResults;

for ne = 1:numel(baseResults)
    cellName = baseResults(ne).name;
    candidateIdx = findCandidateIndices(data, cellName);

    fullWinner = chooseFullFitWinner(data, candidateIdx);
    bestEntry = baseResults(ne);
    bestEntry.frame = validateFrameConsistency(data, candidateIdx, cellName);
    bestEntry = copyFullFitFields(bestEntry, data(fullWinner.dataIdx).results(candidateIdx(fullWinner.dataIdx)));

    nFolds = numel(baseResults(ne).foldLogLikelihood);
    checkCalibration(baseResults(ne), expectedCalibration, data(1).fileName);
    validateFoldCount(baseResults(ne), nFolds, data(1).fileName, cellName);

    for k = 1:nFolds
        foldWinner = chooseFoldWinner(data, candidateIdx, k, nFolds, cellName, expectedCalibration);
        winnerEntry = data(foldWinner.dataIdx).results(candidateIdx(foldWinner.dataIdx));
        bestEntry = copyFoldFields(bestEntry, winnerEntry, k, nFolds);
    end

    bestEntry.cvLogLikelihood = sum(bestEntry.foldLogLikelihood);
    bestEntry.naiveCvLogLikelihood = sum(bestEntry.foldNaiveLogLikelihood);
    bestEntry.cvDeltaLogLikelihood = bestEntry.cvLogLikelihood - bestEntry.naiveCvLogLikelihood;
    bestEntry.cvLogLikelihoodPerFrame = bestEntry.cvLogLikelihood / bestEntry.frame;

    results_best(ne) = bestEntry;
end

for f = 1:numel(groupFiles)
    if ~strcmp(groupFiles(f).calibration, expectedCalibration)
        error('Mixed calibrations found in group %s.', groupFiles(f).group);
    end
end
end

function validateRequiredFields(results, fileName)
requiredFields = {'name', 'frame', 'foldLogLikelihood', 'foldNaiveLogLikelihood', ...
    'foldTrainLogLikelihood', 'foldParams', 'foldAlpha', 'foldOffset', ...
    'foldCalibrationMethod', 'cvLogLikelihood', 'loglikelihood', 'AIC', 'BIC', ...
    'calibration', 'params', 'fullFitParams', 'fullFitAlpha', 'fullFitOffset'};

for i = 1:numel(requiredFields)
    if ~isfield(results, requiredFields{i})
        error('File %s is missing required field %s.', fileName, requiredFields{i});
    end
end
end

function candidateIdx = findCandidateIndices(data, cellName)
candidateIdx = zeros(1, numel(data));

for f = 1:numel(data)
    names = {data(f).results.name};
    idx = find(strcmp(names, cellName), 1);
    if isempty(idx)
        error('Cell %s was not found in %s.', cellName, data(f).fileName);
    end
    candidateIdx(f) = idx;
end
end

function frame = validateFrameConsistency(data, candidateIdx, cellName)
frame = data(1).results(candidateIdx(1)).frame;
if isempty(frame) || ~isscalar(frame) || ~isfinite(frame) || frame <= 0
    error('Invalid frame value for %s in %s.', cellName, data(1).fileName);
end

for f = 2:numel(data)
    thisFrame = data(f).results(candidateIdx(f)).frame;
    if isempty(thisFrame) || ~isscalar(thisFrame) || ~isfinite(thisFrame) || thisFrame <= 0
        error('Invalid frame value for %s in %s.', cellName, data(f).fileName);
    end
    if thisFrame ~= frame
        error('Frame mismatch for %s: %s has %g, %s has %g.', ...
            cellName, data(1).fileName, frame, data(f).fileName, thisFrame);
    end
end
end

function winner = chooseFullFitWinner(data, candidateIdx)
winner.dataIdx = 1;
bestValue = data(1).results(candidateIdx(1)).loglikelihood;

for f = 2:numel(data)
    value = data(f).results(candidateIdx(f)).loglikelihood;
    if value > bestValue
        bestValue = value;
        winner.dataIdx = f;
    end
end
end

function winner = chooseFoldWinner(data, candidateIdx, foldIdx, nFolds, cellName, expectedCalibration)
winner.dataIdx = 1;
bestValue = -Inf;

for f = 1:numel(data)
    entry = data(f).results(candidateIdx(f));
    checkCalibration(entry, expectedCalibration, data(f).fileName);
    validateFoldCount(entry, nFolds, data(f).fileName, cellName);
    value = entry.foldTrainLogLikelihood(foldIdx); % Select fold winner by train LL to avoid data leak.
    if value > bestValue
        bestValue = value;
        winner.dataIdx = f;
    end
end
end

function checkCalibration(entry, expectedCalibration, fileName)
if ~strcmp(char(entry.calibration), expectedCalibration)
    error('Calibration mismatch in %s: expected %s, found %s.', ...
        fileName, expectedCalibration, char(entry.calibration));
end
end

function validateFoldCount(entry, nFolds, fileName, cellName)
foldFields = {'foldLogLikelihood', 'foldNaiveLogLikelihood', 'foldTrainLogLikelihood', ...
    'foldAlpha', 'foldOffset', 'foldCalibrationMethod'};

for i = 1:numel(foldFields)
    if numel(entry.(foldFields{i})) ~= nFolds
        error('Fold count mismatch for %s in %s field %s.', cellName, fileName, foldFields{i});
    end
end

if size(entry.foldParams, 1) ~= nFolds
    error('Fold count mismatch for %s in %s field foldParams.', cellName, fileName);
end
end

function bestEntry = copyFoldFields(bestEntry, winnerEntry, foldIdx, nFolds)
foldFields = {'foldLogLikelihood', 'foldNaiveLogLikelihood', 'foldTrainLogLikelihood', ...
    'foldParams', 'foldAlpha', 'foldOffset', 'foldCalibrationMethod'};

for i = 1:numel(foldFields)
    fieldName = foldFields{i};
    bestEntry.(fieldName) = copyOneFold(bestEntry.(fieldName), winnerEntry.(fieldName), foldIdx, nFolds);
end
end

function dest = copyOneFold(dest, source, foldIdx, nFolds)
if isvector(source)
    dest(foldIdx) = source(foldIdx);
elseif size(source, 1) == nFolds
    dest(foldIdx, :) = source(foldIdx, :);
else
    error('Cannot copy fold %d from field with unsupported size.', foldIdx);
end
end

function bestEntry = copyFullFitFields(bestEntry, winnerEntry)
fullFields = {'loglikelihood', 'AIC', 'BIC', 'params', 'fullFitParams', 'fullFitAlpha', ...
    'fullFitOffset'};

for i = 1:numel(fullFields)
    fieldName = fullFields{i};
    if isfield(bestEntry, fieldName) && isfield(winnerEntry, fieldName)
        bestEntry.(fieldName) = winnerEntry.(fieldName);
    end
end
end
