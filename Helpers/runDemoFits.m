function runDemoFits(demoData, config)
%RUNDEMOFITS Fit every configured RF model and paired calibration.

for runId = config.runIds
    fprintf('\n=== Fit run %d/%d ===\n', runId, config.nRuns);
    for rfIdx = 1:numel(config.rfTypes)
        rfType = config.rfTypes{rfIdx};
        freeFile = resultFile(config, rfType, 'glm_free', runId);
        posFile = resultFile(config, rfType, 'glm_pos', runId);

        if ~config.overwriteFits && isfile(freeFile) && isfile(posFile)
            fprintf('Skipping existing paired fit: %s, run %d\n', ...
                rfType, runId);
            continue
        end

        fprintf('Fitting %s with paired glm_free/glm_pos calibration\n', rfType);
        fitOptions = config.fitOptions;
        fitOptions.UseParallel = config.useParallel;
        fitOptions.SeedTaskId = runId;
        fitOptions.SeedNeuronName = demoData.name;

        [freeResult, posResult] = fitClosestRfCvGlmPair( ...
            demoData.egoClosest, demoData.fired, rfType, ...
            config.numClosest, demoData.folds, ...
            demoData.foldNaiveLogLikelihood, fitOptions);
        freeResult.name = demoData.name;
        freeResult.sparseness = demoData.sparseness;
        posResult.name = demoData.name;
        posResult.sparseness = demoData.sparseness;

        results_to_save = freeResult; %#ok<NASGU>
        save(freeFile, 'results_to_save');
        results_to_save = posResult; %#ok<NASGU>
        save(posFile, 'results_to_save');
    end
end
end

function filePath = resultFile(config, rfType, calibration, runId)
filePath = fullfile(config.resultDir, sprintf( ...
    '%s_%dclosest_cv_%s_%d.mat', ...
    rfType, config.numClosest, calibration, runId));
end
