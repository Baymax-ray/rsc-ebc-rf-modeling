function validateDemoEnvironment(config)
%VALIDATEDEMOENVIRONMENT Check paths, toolboxes, and helper resolution.

if ~isfolder(config.dataDir)
    error('Demo Data folder not found: %s', config.dataDir);
end
if exist('ga', 'file') ~= 2
    error('Global Optimization Toolbox is required (ga was not found).');
end
if exist('fmincon', 'file') ~= 2
    error('Optimization Toolbox is required (fmincon was not found).');
end
if ~ispc
    warning('This demo is designed and tested for Windows MATLAB.');
end

if ~isfolder(config.resultDir)
    mkdir(config.resultDir);
end
if ~isfolder(config.figureDir)
    mkdir(config.figureDir);
end

requiredHelpers = { ...
    'prepareDemoData', ...
    'fitClosestRfCvGlmPair', ...
    'selectDemoBestFits', ...
    'runDemoMonteCarlo', ...
    'drawDemoFigures'};
helperRoot = char(java.io.File(config.helperDir).getCanonicalPath());
helperPrefix = [helperRoot filesep];

for i = 1:numel(requiredHelpers)
    resolved = which(requiredHelpers{i});
    if isempty(resolved)
        error('Required helper was not found: %s', requiredHelpers{i});
    end
    resolved = char(java.io.File(resolved).getCanonicalPath());
    if ~strncmpi(resolved, helperPrefix, numel(helperPrefix))
        error('Helper %s resolved outside the Demo folder: %s', ...
            requiredHelpers{i}, resolved);
    end
end
end
