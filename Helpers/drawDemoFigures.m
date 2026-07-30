function drawDemoFigures(demoData, config)
%DRAWDEMOFIGURES Export component figures and the assembled Figure 2.

basicDir = fullfile(config.figureDir, 'figures_basic');
modelRoot = fullfile(config.figureDir, ...
    sprintf('%dpts Pos', config.numClosest));
if ~isfolder(basicDir)
    mkdir(basicDir);
end
if ~isfolder(modelRoot)
    mkdir(modelRoot);
end

name = demoData.name;
naiveStream = RandStream('mt19937ar', 'Seed', config.naiveSeed);
naiveFired = rand(naiveStream, demoData.frame, 1) < demoData.sparseness;

fig = spatialFigure(demoData, demoData.fired, ...
    sprintf('Recorded Firing (fire rate = %.3f Hz)', ...
    demoData.sparseness * config.fps));
exportDemoFigure(fig, fullfile(basicDir, [name '_raw_firing']));
close(fig);

fig = spatialFigure(demoData, naiveFired, ...
    sprintf('Naive Regenerated Firing (fire rate = %.3f Hz)', ...
    mean(naiveFired) * config.fps));
exportDemoFigure(fig, fullfile( ...
    basicDir, [name '_naive_regenerated_firing']));
close(fig);

allModelSpecs = figureModelSpecs();
configuredTypes = string(config.rfTypes);
modelSpecs = allModelSpecs(ismember(string({allModelSpecs.rfType}), ...
    configuredTypes));
models = repmat(struct(), numel(modelSpecs), 1);
for i = 1:numel(modelSpecs)
    spec = modelSpecs(i);
    filledFile = fullfile(config.resultDir, sprintf( ...
        '%s_%dclosest_cv_%s_filled.mat', ...
        spec.rfType, config.numClosest, config.figureCalibration));
    if ~isfile(filledFile)
        error('Filled result required for figure drawing: %s', filledFile);
    end

    S = load(filledFile, 'results_best');
    if numel(S.results_best) ~= 1
        error('The Demo expects one result in %s.', filledFile);
    end
    entry = S.results_best(1);
    if ~strcmp(char(string(entry.calibration)), config.figureCalibration)
        error('Figure calibration mismatch in %s.', filledFile);
    end

    [RF, xRange, yRange] = computeRfMap(entry, spec.rfType);
    [regenerated, firingProbability] = regeneratedFiring( ...
        demoData, entry, spec.rfType, config.regeneratedSeed);
    rate = mean(regenerated) * config.fps;
    if abs(rate - entry.reFireRate) > 1e-12
        error('Regenerated firing rate does not match the filled result.');
    end

    modelDir = fullfile(modelRoot, spec.rfType);
    if ~isfolder(modelDir)
        mkdir(modelDir);
    end

    fig = rfFigure(RF, xRange, yRange, entry.fullFitOffset);
    rfBase = sprintf('%s_%s_%dclosest_rf', ...
        name, spec.rfType, config.numClosest);
    exportDemoFigure(fig, fullfile(modelDir, rfBase));
    close(fig);

    fig = spatialFigure(demoData, regenerated, sprintf( ...
        'Regenerated Firing (fire rate = %.3f Hz, p = %.4g)', ...
        rate, entry.pvalue));
    firingBase = sprintf('%s_%s_%dclosest_regenerated_firing', ...
        name, spec.rfType, config.numClosest);
    exportDemoFigure(fig, fullfile(modelDir, firingBase));
    close(fig);

    models(i).label = spec.label;
    models(i).title = spec.title;
    models(i).rfType = spec.rfType;
    models(i).entry = entry;
    models(i).RF = RF;
    models(i).xRange = xRange;
    models(i).yRange = yRange;
    models(i).regenerated = regenerated;
    models(i).firingProbability = firingProbability;
end

fig = figure('Visible', 'off', 'Color', 'white', ...
    'Units', 'centimeters', 'Position', [2, 2, 6.5, 6.5]);
ax = axes(fig, 'Position', [0.08, 0.08, 0.84, 0.80]);
drawHeadDirectionCompass(ax, true);
exportDemoFigure(fig, fullfile(config.figureDir, 'directional_compass'));
close(fig);

if numel(models) == numel(allModelSpecs)
    drawDemoFigure2(demoData, naiveFired, models, config);
else
    fprintf('Skipping assembled Figure 2 because only %d of %d RF types are configured.\n', ...
        numel(models), numel(allModelSpecs));
end
end

function fig = spatialFigure(demoData, fired, titleText)
fig = figure('Visible', 'off', 'Color', 'white', ...
    'Units', 'centimeters', 'Position', [2, 2, 10, 9]);
ax = axes(fig);
scatter(ax, demoData.xPosition, demoData.yPosition, ...
    1 + 20 * double(fired), demoData.headDirection, 'filled');
colormap(ax, hsv(360));
caxis(ax, [0, 360]);
xlim(ax, [0, 120]);
ylim(ax, [0, 120]);
axis(ax, 'square');
xlabel(ax, 'X position (cm)');
ylabel(ax, 'Y position (cm)');
title(ax, titleText);
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off');
end

function fig = rfFigure(RF, xRange, yRange, offset)
fig = figure('Visible', 'off', 'Color', 'white', ...
    'Units', 'centimeters', 'Position', [2, 2, 10, 9]);
ax = axes(fig);
imagesc(ax, xRange, yRange, RF);
axis(ax, 'xy');
axis(ax, 'square');
applyRfLogitContributionStyle(fig, RF, offset);
xlabel(ax, 'Right (+) / Left (-)  [cm]');
ylabel(ax, 'Forward (+) / Backward (-)  [cm]');
hold(ax, 'on');
plot(ax, 0, 0, 'ko', 'MarkerSize', 5, 'MarkerFaceColor', 'black');
quiver(ax, 0, 0, 0, 10, 'black', 'LineWidth', 1.5, ...
    'MaxHeadSize', 1.5);
hold(ax, 'off');
set(ax, 'FontName', 'Arial', 'FontSize', 10, 'Box', 'off');
end

function [RF, xRange, yRange] = computeRfMap(entry, rfType)
xRange = linspace(-30, 30, 500);
yRange = linspace(-30, 30, 500);
[X, Y] = meshgrid(xRange, yRange);
egoGrid = zeros(numel(X), 1, 2);
egoGrid(:, 1, 1) = Y(:);
egoGrid(:, 1, 2) = -X(:);
rfValue = computeRfValueForClosest( ...
    rfType, egoGrid, entry.fullFitParams);
RF = reshape(entry.fullFitAlpha * rfValue, size(X));
end

function [regenerated, probability] = regeneratedFiring( ...
    demoData, entry, rfType, seed)
rfValue = computeRfValueForClosest( ...
    rfType, demoData.egoClosest, entry.fullFitParams);
eta = entry.fullFitAlpha * rfValue - entry.fullFitOffset;
probability = 1 ./ (1 + exp(-eta));
epsValue = 1e-15;
probability = min(max(probability, epsValue), 1 - epsValue);
stream = RandStream('mt19937ar', 'Seed', seed);
regenerated = rand(stream, numel(probability), 1) < probability;
end

function specs = figureModelSpecs()
specs = struct( ...
    'label', {'C', 'D', 'E', 'F', 'G', 'H', 'I'}, ...
    'title', { ...
        'Angle-Distance Gaussian', ...
        'Angle-only Gaussian', ...
        'Distance-only Gaussian', ...
        'Isotropic 2D Gaussian', ...
        'Elliptical 2D Gaussian', ...
        'Angle-bounded Distance Gaussian', ...
        'Difference of Gaussians (DoG)'}, ...
    'rfType', { ...
        'angle_distance_gaussian', ...
        'angle_gaussian', ...
        'distance_gaussian', ...
        'single_2d_gaussian', ...
        'stretched_2d_gaussian', ...
        'bounded_gaussian', ...
        'dog'});
end
