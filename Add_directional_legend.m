function Add_directional_legend(ax_position)
% Create the compass like directional legend for a plot
%
% Inputs:
%   ax_position: 1x4 vector specifying the position of the legend axes ([x, y, width, height])

    % Create small axes
    legend_ax = axes('Position', ax_position);

    % Generate annular (ring-shaped) color wheel — exclude center
    theta = linspace(0, 2*pi, 361);
    r = linspace(0.5, 1, 2);  % avoid r = 0 to leave center blank
    [THETA,R] = meshgrid(theta, r);
    X = R .* cos(THETA);
    Y = R .* sin(THETA);
    C = rad2deg(THETA);  % color by angle in degrees

    % Plot using pcolor
    p = pcolor(legend_ax, X, Y, C);
    p.EdgeColor = 'none';
    shading(legend_ax, 'flat');
    colormap(legend_ax, hsv);
    clim(legend_ax, [0 360]);

    axis(legend_ax, 'equal', 'off');

    % Add directional labels
    text(legend_ax, 1.3, 0, '0°',   'HorizontalAlignment','center', 'FontSize',8, 'FontWeight','bold');
    text(legend_ax, 0, 1.3, '90°',  'HorizontalAlignment','center', 'FontSize',8, 'FontWeight','bold');
    text(legend_ax, -1.6, 0, '180°','HorizontalAlignment','center', 'FontSize',8, 'FontWeight','bold');
    text(legend_ax, 0, -1.2, '270°','HorizontalAlignment','center', 'FontSize',8, 'FontWeight','bold');


    text(legend_ax, 0, 0, 'HD', 'HorizontalAlignment','center', 'FontSize',7, 'FontWeight','bold');
end