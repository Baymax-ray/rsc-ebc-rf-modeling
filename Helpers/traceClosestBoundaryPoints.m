function [Egorelative_closest, dist_all, intersectionPts_all] = ...
    traceClosestBoundaryPoints(x_position, y_position, head_direction, ...
    boundaries_all, numClosest, angles, useParallel)
%TRACECLOSESTBOUNDARYPOINTS Trace boundary hits with a serial fallback.

if nargin < 6 || isempty(angles)
    angles = 0:3:357;
end
if nargin < 7 || isempty(useParallel)
    useParallel = false;
end
if numClosest > numel(angles)
    error('numClosest cannot exceed the number of sampled ray angles.');
end

x_position = x_position(:);
y_position = y_position(:);
head_direction = head_direction(:);
nt = numel(x_position);

Egorelative_closest = zeros(nt, numClosest, 2);
dist_all = zeros(nt, numel(angles));
intersectionPts_all = zeros(nt, numel(angles), 2);

if useParallel
    anglesConst = parallel.pool.Constant(angles);
    boundariesConst = parallel.pool.Constant(boundaries_all);
    parfor t = 1:nt
        [ego, dist, intersections] = traceOneFrame( ...
            x_position(t), y_position(t), head_direction(t), ...
            boundariesConst.Value, numClosest, anglesConst.Value, t);
        Egorelative_closest(t, :, :) = ego;
        dist_all(t, :) = dist;
        intersectionPts_all(t, :, :) = intersections;
    end
else
    for t = 1:nt
        [ego, dist, intersections] = traceOneFrame( ...
            x_position(t), y_position(t), head_direction(t), ...
            boundaries_all, numClosest, angles, t);
        Egorelative_closest(t, :, :) = ego;
        dist_all(t, :) = dist;
        intersectionPts_all(t, :, :) = intersections;
    end
end
end

function [ego, dist, intersectionPts] = traceOneFrame( ...
    x, y, headDirection, boundaries_all, numClosest, angles, frameIdx)
position = [x, y];
boundaries = Update_boundaries(boundaries_all, frameIdx);
[dist, intersectionPts] = multiRayTrace2D( ...
    position, headDirection + angles, boundaries);

[~, order] = sort(dist);
closestPoints = intersectionPts(order(1:numClosest), :);
R = [cosd(headDirection), sind(headDirection); ...
    -sind(headDirection), cosd(headDirection)];
ego = (R * (closestPoints - position)')';
end
