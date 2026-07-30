function [distances, intersectionPts] = multiRayTrace2D(P, angles, walls)
% MULTIRAYTRACE2D  Vectorized version of Ray_Trace2D for multiple rays.
%
%   P         : [px, py], the starting point of the rays.
%   angles    : Mx1 array of angles in degrees.
%    walls    : struct array with fields:
%       .type ('segment' or 'arc')
%       .x1, .y1, .x2, .y2 (for segments)
%       .cx, .cy, .radius, .start_angle, .end_angle (for arcs)
%       .appear_timestamp, .disappear_timestamp (not needed here)
%   distances : Mx1 array of the *closest* intersection distance for each angle.
%               If no intersection, the distance is Inf.
%   intersectionPt: Mx2 array of indices indicating the intersection point [x,y] on the wall segment. If none is hit, that entry is [Inf, Inf].
%
% Example:
%   walls = struct('type', 'segment', 'x1', 0, 'y1', 0, 'x2', 10, 'y2', 0);
%   % 72 angles in 5-degree increments
%   angles = 0:5:355;
%   P = [60, 40];
%   [dist, wIdx] = multiRayTrace2D(P, angles, walls);

    % Ensure angles is a column vector for consistent broadcasting
    angles = angles(:);  % Mx1

    M = numel(angles);   % number of rays
    N = size(walls,1);   % number of walls

    % Precompute direction vectors [dX, dY] for each angle (size Mx2).
    dX = cosd(angles);   % Mx1
    dY = sind(angles);   % Mx1
    
    % Unit direction vectors (already unit if using cosd/sind)
    % but let's be explicit if you like:
    % norms = sqrt(dX.^2 + dY.^2);
    % dX = dX ./ norms;
    % dY = dY ./ norms;

    % Ray start
    px = P(1);
    py = P(2);
    
    x1 = zeros(N,1);
    y1 = zeros(N,1);
    x2 = zeros(N,1);
    y2 = zeros(N,1);
    % Extract wall endpoints (size Nx1 for each coordinate)
    for i = 1:length(walls)
        w= walls(i);
        switch w.type
            case "segment"
                x1(i) = w.x1;
                y1(i) = w.y1;
                x2(i) = w.x2;
                y2(i) = w.y2;
            otherwise
                error('Unknown wall type: %s', w.type);
        end
    end
    % Wall direction vectors
    Wx = x2 - x1;  % Nx1
    Wy = y2 - y1;  % Nx1

    % We will compute a 2D "grid" of intersections:
    %   each row corresponds to one angle (ray),
    %   each column corresponds to one wall.
    % Denominator for each (ray, wall): size MxN
    % denom = (dX * Wy') - (dY * Wx')
    % We'll do it in a way that results in MxN:
    denom = dX .* Wy' - dY .* Wx';

    % Differences from the ray start to each wall's start:
    dx = (x1 - px);  % Nx1
    dy = (y1 - py);  % Nx1

    % We also need these in MxN, so do the same broadcast:
    %    dx' is 1xN, dX is Mx1 => (M x N)
    % t = [ (dx * Wy - dy * Wx) / denom ]
    numerT = dx .* Wy - dy .* Wx;  % Nx1
    % Now replicate across M rows => (M x N)
    numerT = repmat(numerT', M, 1);

    t = numerT ./ denom;  % MxN

    % Similarly for u:
    % u = [((x1-px)*dY - (y1-py)*dX) / denom]
    % But we have M different dY. So let's do the broadcast properly:
    % We want MxN => best approach is to make:
    %   dx:  Nx1   => replicate to MxN
    %   dY:  Mx1   => replicate to MxN
    % We'll do:
    dx2 = repmat(dx', M, 1);   % MxN
    dy2 = repmat(dy', M, 1);   % MxN
    dX2 = repmat(dX, 1, N);    % MxN
    dY2 = repmat(dY, 1, N);    % MxN

    numerU = (dx2 .* dY2) - (dy2 .* dX2);  % MxN

    % Now t = numerT ./ denom (already MxN),
    % and u = numerU ./ denom (MxN).
    u = numerU ./ denom;

    t= round(t,12); %to avoid floating point error
    u= round(u,12); %to avoid floating point error
    % Mark parallel or nearly-parallel as invalid => denom ~ 0
    parallelMask = abs(denom) < 1e-12;
    t(parallelMask) = Inf;
    u(parallelMask) = Inf;

    % Valid intersection => t >= 0,  0 <= u <= 1
    % We'll create a logical mask for each (angle,wall).
    validMask = (t >= 0) & (u >= 0) & (u <= 1);

    % If an entry is invalid, set t to Inf so it won't be chosen as a minimum
    t(~validMask) = Inf;

    % Now for each row (angle), we find the minimal distance across all walls
    minT= min(t, [], 2);

    distances = minT;      % Mx1
    if any(isinf(distances))
        disp('Warning: some rays did not intersect any walls.');
    end
    intersectionPts= [px + minT .* dX, py + minT .* dY];% Mx1
    intersectionPts(isnan(intersectionPts)) = Inf; %NaN is created by Inf-Inf
end
