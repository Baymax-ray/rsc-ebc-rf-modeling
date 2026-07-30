%% Patrick nearest-boundary CV/GLM demo
% This self-contained Windows MATLAB demo reproduces the complete analysis
% sequence for the example neuron used in Figure 2:
%   data preparation -> repeated fitting -> best selection ->
%   Monte Carlo validation -> component and assembled figure export.
%
% The Demo does not use SCC shell scripts, SGE_TASK_ID, environment
% variables, or code outside this folder.

%% 1. Configuration and Windows setup
% Use "paper" for the full analysis (5 runs, 5 folds, 2000 null draws).
% Use "quick" for a shorter workflow demonstration that is not intended to
% reproduce the final paper values.
mode = "paper";

demoRoot = fileparts(mfilename('fullpath'));
if isempty(demoRoot)
    error('Run Demo.m from a saved copy of the code and demo folder.');
end
helperDir = fullfile(demoRoot, 'Helpers');
addpath(helperDir, '-begin');
rehash path;

config = makeDemoConfig(demoRoot, mode);

% Stage switches allow expensive fits to be reused. When fit is false, the
% configured numeric run files must already exist in Demo Results/<mode>.
config.stages.prepareData = true;
config.stages.fit = true;
config.stages.selectBest = true;
config.stages.monteCarlo = true;
config.stages.drawFigures = true;

% Set overwriteFits=false to reuse complete paired free/positive run files.
config.overwriteFits = true;
config.nWorkers = 4;

validateDemoEnvironment(config);
config.useParallel = configureDemoParallelPool(config.nWorkers);
config.fitOptions.UseParallel = config.useParallel;

fprintf('\nPatrick Demo mode: %s\n', config.mode);
fprintf('Results: %s\n', config.resultDir);
fprintf('Figures: %s\n', config.figureDir);
save(config.manifestFile, 'config');

%% 2. Load, validate, and prepare the example neuron
% Spikes are assigned to the nearest tracking frame and represented as a
% Bernoulli response. Closest visible boundary points are ray-traced once
% and reused by every model, run, fold, and Monte Carlo calculation.
if config.stages.prepareData
    demoData = prepareDemoData(config);
    save(config.preparedDataFile, 'demoData', '-v7.3');
else
    if ~isfile(config.preparedDataFile)
        error('Prepared data cache not found: %s', config.preparedDataFile);
    end
    cached = load(config.preparedDataFile, 'demoData');
    demoData = cached.demoData;
end
validatePreparedDemoData(demoData, config);

%% 3. Repeated blocked-CV fitting
% RF shape parameters are optimized only on each training partition.
% For a fixed shape, alpha and offset are fitted by logistic calibration on
% the same training frames. Held-out frames are used only for CV scoring.
% glm_free and glm_pos are fitted as a pair and saved separately.
if config.stages.fit
    runDemoFits(demoData, config);
end

%% 4. Per-fold and full-fit best selection
% Fold winners are selected by foldTrainLogLikelihood to avoid data leak.
% Held-out fold likelihoods from those winners are stitched into the final
% CV score. Full-fit fields, including AIC and BIC, independently follow
% the run with the largest full-data loglikelihood.
if config.stages.selectBest
    selectDemoBestFits(config);
end

%% 5. Monte Carlo CV validation and regenerated firing
% Stored fold parameters reconstruct fixed out-of-fold probabilities.
% Fold-specific mean-rate Bernoulli data form the null distribution.
% Full-fit probabilities generate a deterministic regenerated firing map.
if config.stages.monteCarlo
    runDemoMonteCarlo(demoData, config);
end

%% 6. Component figures and assembled Figure 2
% Individual RF maps display alpha*g(theta) without offset or sigmoid.
% The logistic zero-input baseline is reported in each RF subtitle.
% Figure 2 uses glm_pos results and the manually adjusted A-I layout.
if config.stages.drawFigures
    drawDemoFigures(demoData, config);
end

%% 7. Completion summary
fprintf('\nDemo complete.\n');
fprintf('MAT results: %s\n', config.resultDir);
fprintf('PDF/PNG figures: %s\n', config.figureDir);
