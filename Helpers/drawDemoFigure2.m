function drawDemoFigure2(demoData, naiveFired, models, config)
%DRAWDEMOFIGURE2 Assemble the self-contained MATLAB version of Figure 2.

pageWidth = 170;
pageHeight = 235;
fig = figure('Visible', 'off', 'Color', 'white', ...
    'Units', 'centimeters', 'Position', [1, 1, 17, 23.5], ...
    'PaperUnits', 'centimeters', 'PaperSize', [17, 23.5], ...
    'PaperPosition', [0, 0, 17, 23.5]);

drawTopPanel(fig, demoData, demoData.fired, ...
    'A', 'Recorded firing', 4, true, pageWidth, pageHeight);
drawTopPanel(fig, demoData, naiveFired, ...
    'B', 'Naive baseline', 87, false, pageWidth, pageHeight);

rows = [0, 0, 1, 1, 2, 2, 3];
columns = [0, 1, 0, 1, 0, 1, 0];
for i = 1:numel(models)
    drawModelBlock(fig, demoData, models(i), rows(i), columns(i), ...
        i == 1, pageWidth, pageHeight);
end
drawCompassBlock(fig, pageWidth, pageHeight);

basePath = fullfile(config.figureDir, 'figure2_N165_nearest_only');
exportDemoFigure(fig, basePath, true);
close(fig);
end

function drawTopPanel(fig, demoData, fired, label, titleText, ...
    blockX, showAxes, pageWidth, pageHeight)
blockWidth = 79;
plotSide = 34.5;
plotX = blockX + (blockWidth - plotSide) / 2;
plotY = 11.5;

addHeading(fig, label, titleText, blockX, 4.5, pageWidth, pageHeight, 9.5);
ax = axes(fig, 'Units', 'normalized', ...
    'Position', topRect(plotX, plotY, plotSide, plotSide, ...
    pageWidth, pageHeight));
plotSpatialMap(ax, demoData, fired);
if showAxes
    set(ax, 'XTick', [0, 60, 120], 'YTick', [0, 60, 120], ...
        'FontSize', 8);
    xlabel(ax, 'X position (cm)', 'FontName', 'Arial', 'FontSize', 9);
    ylabel(ax, 'Y position (cm)', 'FontName', 'Arial', 'FontSize', 9);
else
    set(ax, 'XTick', [], 'YTick', []);
end
end

function drawModelBlock(fig, demoData, model, row, column, highlight, ...
    pageWidth, pageHeight)
if column == 0
    blockX = 4;
else
    blockX = 87;
end
blockY = 55 + row * 45;
blockWidth = 79;
blockHeight = 44;

if highlight
    axes(fig, 'Units', 'normalized', ...
        'Position', topRect(blockX, blockY, blockWidth, blockHeight, ...
        pageWidth, pageHeight), ...
        'Color', [0.918, 0.957, 0.984], ...
        'XColor', [0.216, 0.459, 0.729], ...
        'YColor', [0.216, 0.459, 0.729], ...
        'LineWidth', 1, 'Box', 'on', ...
        'XTick', [], 'YTick', [], 'HitTest', 'off');
end

addHeading(fig, model.label, model.title, blockX, blockY + 2.5, ...
    pageWidth, pageHeight, 8.8);

plotSide = 24.1;
plotY = blockY + 13;
if highlight
    rfX = blockX + 10;
else
    rfX = blockX + 4.5;
end
spatialX = blockX + 47;
colorbarX = rfX + plotSide + 0.9;

addCenteredText(fig, 'RF contribution', rfX, blockY + 9.3, plotSide, ...
    pageWidth, pageHeight, 8, [0.43, 0.43, 0.43]);
addCenteredText(fig, 'Regenerated firing', spatialX, blockY + 9.3, ...
    plotSide, pageWidth, pageHeight, 8, [0.43, 0.43, 0.43]);

rfPosition = topRect(rfX, plotY, plotSide, plotSide, ...
    pageWidth, pageHeight);
axRf = axes(fig, 'Units', 'normalized', 'Position', rfPosition);
imagesc(axRf, model.xRange, model.yRange, model.RF);
axis(axRf, 'xy');
axis(axRf, 'square');
set(axRf, 'Position', rfPosition, 'FontName', 'Arial', ...
    'FontSize', 7, 'Box', 'off', 'Layer', 'top');
[colorLimits, colorMap, isConstant] = ...
    getSignedRfColorSettings(model.RF);
colormap(axRf, colorMap);
caxis(axRf, colorLimits);
hold(axRf, 'on');
plot(axRf, 0, 0, 'ko', 'MarkerSize', 3, 'MarkerFaceColor', 'black');
quiver(axRf, 0, 0, 0, 10, 0, 'Color', 'black', ...
    'LineWidth', 0.8, 'MaxHeadSize', 1.5);
hold(axRf, 'off');

if highlight
    set(axRf, 'XTick', [-30, 0, 30], 'YTick', [-30, 0, 30]);
    xlabel(axRf, 'Left (-)       cm       Right (+)', ...
        'FontName', 'Arial', 'FontSize', 7);
    ylabel(axRf, 'Backward (-) / Forward (+)', ...
        'FontName', 'Arial', 'FontSize', 7);
else
    set(axRf, 'XTick', [], 'YTick', []);
end

cb = colorbar(axRf);
set(cb, 'Units', 'normalized', ...
    'Position', topRect(colorbarX, plotY, 1.15, plotSide, ...
    pageWidth, pageHeight), ...
    'FontName', 'Arial', 'FontSize', 7);
if isConstant
    cb.Ticks = model.RF(1);
else
    cb.Ticks = linspace(min(model.RF(:)), max(model.RF(:)), 3);
end
cb.TickLabels = arrayfun(@(value) sprintf('%.3f', value), ...
    cb.Ticks, 'UniformOutput', false);
set(axRf, 'Position', rfPosition);

spatialPosition = topRect(spatialX, plotY, plotSide, plotSide, ...
    pageWidth, pageHeight);
axSpatial = axes(fig, 'Units', 'normalized', 'Position', spatialPosition);
plotSpatialMap(axSpatial, demoData, model.regenerated);
set(axSpatial, 'XTick', [], 'YTick', [], 'Position', spatialPosition);
end

function drawCompassBlock(fig, pageWidth, pageHeight)
blockX = 87;
blockY = 190;
blockWidth = 79;
addCenteredText(fig, 'Spike colour denotes head direction', ...
    blockX, blockY + 1.8, blockWidth, pageWidth, pageHeight, ...
    9.5, [0.145, 0.145, 0.145], 'bold');
ax = axes(fig, 'Units', 'normalized', ...
    'Position', topRect(blockX + 19, blockY + 10, 41, 34, ...
    pageWidth, pageHeight));
drawHeadDirectionCompass(ax, false, 200);
end

function plotSpatialMap(ax, demoData, fired)
scatter(ax, demoData.xPosition, demoData.yPosition, ...
    0.2 + 3 * double(fired), demoData.headDirection, 'filled');
colormap(ax, hsv(360));
caxis(ax, [0, 360]);
xlim(ax, [0, 120]);
ylim(ax, [0, 120]);
axis(ax, 'square');
set(ax, 'FontName', 'Arial', 'Box', 'off', 'Layer', 'top');
end

function addHeading(fig, label, titleText, x, y, pageWidth, pageHeight, ...
    titleSize)
annotation(fig, 'textbox', topRect(x, y, 7, 6, pageWidth, pageHeight), ...
    'String', label, 'LineStyle', 'none', 'Margin', 0, ...
    'FontName', 'Arial', 'FontSize', 15, 'Color', [0.145, 0.145, 0.145], ...
    'VerticalAlignment', 'middle');
annotation(fig, 'textbox', topRect(x + 7, y + 0.5, 69, 5, ...
    pageWidth, pageHeight), ...
    'String', titleText, 'LineStyle', 'none', 'Margin', 0, ...
    'FontName', 'Arial', 'FontWeight', 'bold', ...
    'FontSize', titleSize, 'Color', [0.145, 0.145, 0.145], ...
    'VerticalAlignment', 'middle');
end

function addCenteredText(fig, textValue, x, y, width, ...
    pageWidth, pageHeight, fontSize, color, fontWeight)
if nargin < 10
    fontWeight = 'normal';
end
annotation(fig, 'textbox', topRect(x, y, width, 4, ...
    pageWidth, pageHeight), ...
    'String', textValue, 'LineStyle', 'none', 'Margin', 0, ...
    'HorizontalAlignment', 'center', 'VerticalAlignment', 'middle', ...
    'FontName', 'Arial', 'FontWeight', fontWeight, ...
    'FontSize', fontSize, 'Color', color);
end

function position = topRect(x, y, width, height, pageWidth, pageHeight)
position = [x / pageWidth, ...
    1 - (y + height) / pageHeight, ...
    width / pageWidth, ...
    height / pageHeight];
end
