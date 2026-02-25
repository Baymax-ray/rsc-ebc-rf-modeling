function [dist,intersectionPt] = Ray_Trace2D(P, theta, walls)
% RAYTRACE2D  Finds the distance from point P to the closest wall
% in allocentric direction theta (in degrees) among a list of line segments ('walls').
%
%   P        : [px, py], the starting point of the ray
%   theta    : angle in degrees measured from the x-axis
%   walls    : struct array with fields:
%       .type ('segment' or 'arc')
%       .x1, .y1, .x2, .y2 (for segments)
%       .cx, .cy, .radius, .start_angle, .end_angle (for arcs)
%       .appear_timestamp, .disappear_timestamp (not needed here)
%   dist     : distance to the nearest intersection. If no wall is hit, returns Inf.
%   intersectionPt : [ix, iy], the intersection point on the wall segment
%
% Example usage:
%   walls = struct('type', 'segment', 'x1', 0, 'y1', 0, 'x2', 10, 'y2', 0);
%   [dist, intersectionPt] = Ray_Trace2D([10,10], 45, walls);

    % Convert angle to a direction vector
    dX = cosd(theta);
    dY = sind(theta);

    D = [dX; dY];
    normD = norm(D); %actually normD should be 1, just for checking
    if normD < 1e-12
        dist = Inf;
        intersectionPt = [Inf, Inf];
        disp('Warning: zero-length direction vector');
        return;
    end
    % Keep unit direction to have distance = t directly
    D = D / normD;

    px = P(1);
    py = P(2);

    tMin = Inf;   % We'll track the minimum t over all walls
    intersectionPt = [Inf, Inf];

    % Loop over each wall
    for i = 1:length(walls)
        w= walls(i);
        switch w.type
            case "segment"
            x1 = w.x1;
            y1 = w.y1;
            x2 = w.x2;
            y2 = w.y2;

            % Solve for intersection of:
            %   P + t*D = W1 + u*(W2 - W1)
            % in param form:
            %   px + t*dX = x1 + u*(x2 - x1)
            %   py + t*dY = y1 + u*(y2 - y1)

            Wx = x2 - x1;
            Wy = y2 - y1;

            denom = (D(1)*Wy - D(2)*Wx);

            % If denom = 0, lines are parallel or coincident => no single intersection
            if abs(denom) < 1e-12
                continue;
            end

            % Numerators for t and u:
            % Using determinant / cross-product style
            dx = x1 - px;
            dy = y1 - py;

            t = (dx * Wy - dy * Wx) / denom;   % for ray
            u = (dx * D(2) - dy * D(1)) / denom;  % for wall

            % We want t >= 0 (forward along ray)
            % and u between 0 and 1 (within wall segment)
            u = round(u,12);
            t = round(t,12); %to avoid floating point error
            if (t >= 0) && (u >= 0) && (u <= 1)
                % Intersection is valid. Check if it's the nearest so far.
                if t < tMin
                    tMin = t;
                    intersectionPt = [px + t*D(1), py + t*D(2)];
                end
            end
            case "arc"
            % currently not implemented
            otherwise
                error('Unknown wall type: %s', w.type);
        end
    end
    % If tMin is Inf, there was no intersection
    if isinf(tMin)
        dist = Inf;
        disp('Warning: no intersection found');
    else
        % Because D was normalized, tMin is the direct distance
        dist = tMin;
    end

end
