function config = makeDemoConfig(demoRoot, mode)
%MAKEDEMOCONFIG Build the reproducible configuration for the Patrick demo.

if nargin < 2 || isempty(mode)
    mode = 'paper';
end

mode = lower(char(string(mode)));
config = struct();
config.mode = mode;
config.demoRoot = demoRoot;
config.helperDir = fullfile(demoRoot, 'Helpers');
config.dataDir = fullfile(demoRoot, 'Demo Data');
config.resultDir = fullfile(demoRoot, 'Demo Results', mode);
config.figureDir = fullfile(demoRoot, 'Demo Figures', mode);
config.preparedDataFile = fullfile(config.resultDir, 'demo_prepared_data.mat');
config.manifestFile = fullfile(config.resultDir, 'demo_run_manifest.mat');

config.rfTypes = { ...
    'angle_distance_gaussian', ...
    'angle_gaussian', ...
    'distance_gaussian', ...
    'single_2d_gaussian', ...
    'stretched_2d_gaussian', ...
    'dog', ...
    'bounded_gaussian'};
config.calibrations = {'glm_free', 'glm_pos'};
config.figureCalibration = 'glm_pos';
config.numClosest = 1;
config.nWorkers = 4;
config.fps = 30;
config.rayAngles = 0:3:357;

config.baseSeed = 12345;
config.paperNeuronIndex = 36;
config.regeneratedSeed = config.baseSeed + config.paperNeuronIndex;
config.monteCarloSeedBase = ...
    config.baseSeed + config.paperNeuronIndex * 10000;
config.naiveSeed = makeDeterministicSeed( ...
    config.baseSeed, 'N36', 'naive', 'full');

config.stages = struct( ...
    'prepareData', true, ...
    'fit', true, ...
    'selectBest', true, ...
    'monteCarlo', true, ...
    'drawFigures', true);
config.overwriteFits = true;

switch mode
    case 'paper'
        config.nRuns = 5;
        config.nFolds = 5;
        config.nPermutations = 2000;
        config.fitOptions = struct('DoFullFit', true);
    case 'quick'
        config.nRuns = 2;
        config.nFolds = 3;
        config.nPermutations = 200;
        config.fitOptions = struct( ...
            'DoFullFit', true, ...
            'PopulationSize', 80, ...
            'EliteCount', 4, ...
            'MaxGenerations', 20, ...
            'MaxStallGenerations', 8, ...
            'HybridMaxIterations', 200, ...
            'HybridMaxFunctionEvaluations', 3000, ...
            'FunctionTolerance', 1e-4);
    otherwise
        error('Unknown demo mode: %s. Use paper or quick.', mode);
end

config.runIds = 1:config.nRuns;
config.useParallel = false;
config.fitOptions.UseParallel = false;
end
