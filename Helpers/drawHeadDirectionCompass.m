function drawHeadDirectionCompass(ax, showTitle, ringMarkerSize)
%DRAWHEADDIRECTIONCOMPASS Draw the shared HSV head-direction key.

if nargin < 2
    showTitle = true;
end
if nargin < 3
    ringMarkerSize = 14;
end

angles = (0:359)';
x = cosd(angles);
y = sind(angles);
scatter(ax, x, y, ringMarkerSize, angles, 'filled');
colormap(ax, hsv(360));
caxis(ax, [0, 360]);
axis(ax, 'equal');
axis(ax, [-1.35, 1.35, -1.35, 1.35]);
axis(ax, 'off');
hold(ax, 'on');
text(ax, 0, 0, 'HD', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'middle', 'FontName', 'Arial', ...
    'FontWeight', 'bold', 'FontSize', 10);
text(ax, 1.18, 0, '0°', 'HorizontalAlignment', 'left', ...
    'FontName', 'Arial', 'FontSize', 9);
text(ax, 0, 1.17, '90°', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'bottom', 'FontName', 'Arial', 'FontSize', 9);
text(ax, -1.18, 0, '180°', 'HorizontalAlignment', 'right', ...
    'FontName', 'Arial', 'FontSize', 9);
text(ax, 0, -1.17, '270°', 'HorizontalAlignment', 'center', ...
    'VerticalAlignment', 'top', 'FontName', 'Arial', 'FontSize', 9);
if showTitle
    title(ax, 'Spike colour denotes head direction', ...
        'FontName', 'Arial', 'FontWeight', 'bold', 'FontSize', 10);
end
hold(ax, 'off');
end
