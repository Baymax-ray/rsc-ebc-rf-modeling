function val = Check_receptive_field(rftype, position, head_direction, targets, params)
% CHECK_RECEPTIVE_FIELD Return the response of RF using one of these candidate receptive fields for multiple targets at a time.
%
%   val = Check_receptive_field(rftype, position, head_direction, target, params)
%
% Inputs:
%   rftype        : string specifying which RF to use
%                   (e.g., 'angle_distance_gaussian', 'single_2d_gaussian', 'bounded_gaussian')
%   position      : [px, py], the animal's position in 2D
%   head_direction: the animal's heading angle in degrees
%   targets        : N*2 ([tx, ty]), the point at which we evaluate the RF
%   params        : struct with fields depending on rftype
%
% Output:
%   val           : N*scalar, the RF response or weight at the 'target' location.
%
% Example usage:
%   params.sigmaAngle = 30;   % degrees
%   params.sigmaDist  = 5;    % distance units
%   val = Check_receptive_field('angle_distance_gaussian', [0,0], 0, [3,4], params);
%   disp(val);   % e.g. 0.57...

val = 0;
for i=1:size(targets,1)

% -- 1) Compute vector from 'position' to 'target'
dx=targets(i,1)-position(1);
dy=targets(i,2)-position(2);
dist = sqrt(dx^2 + dy^2);        % distance from 'position' to 'target'

% -- 2) Compute angle of that vector, and difference from 'head_direction'
% atan2d returns angle in degrees, range [-180, 180].
angleToTarget = atan2d(dy, dx);
angleDiff = angleToTarget - head_direction;

% Wrap angleDiff into [-180, 180]
angleDiff = mod(angleDiff + 180, 360) - 180;

% -- 3) Evaluate the receptive field depending on rftype
switch lower(rftype)

    case 'angle_distance_gaussian'
        % A 2D Gaussian in "angle × distance" space:
        %   - One Gaussian for the angle difference
        %   - One Gaussian for the radial distance
        % The final RF is the product of these two Gaussians.
        %
        % Required params:
        %   params.sigmaAngle (width in degrees)
        %   params.muAngle    (preferred angle in degrees)
        %   params.sigmaDist  (width in distance units)
        %   params.muDist     (preferred distance in distance units)
        
        sigmaAngle = params.sigmaAngle;
        sigmaDist  = params.sigmaDist;
        angleDeviation = angleDiff - params.muAngle;
        angleDeviation = mod(angleDeviation + 180, 360) - 180;  % wrap to [-180, 180]


        gaussAngle = exp(-0.5 * (angleDeviation^2) / (sigmaAngle^2));
        gaussDist  = exp(-0.5 * ((dist - params.muDist)^2) / (sigmaDist^2));

        val = val+gaussAngle * gaussDist;

    case 'angle_gaussian'
        % A 1D Gaussian in "angle" space:
        %   - One Gaussian for the angle difference
        %
        % Required params:
        %   params.sigmaAngle (width in degrees)
        %   params.muAngle    (preferred angle in degrees)

        sigmaAngle = params.sigmaAngle;
        angleDeviation = angleDiff - params.muAngle;
        angleDeviation = mod(angleDeviation + 180, 360) - 180;  % wrap to [-180, 180]

        gaussAngle = exp(-0.5 * (angleDeviation^2) / (sigmaAngle^2));
        
        val = val+gaussAngle;

    case 'distance_gaussian'
        % A 1D Gaussian in "distance" space:
        %   - One Gaussian for the radial distance
        %
        % Required params:
        %   params.sigmaDist  (width in distance units)
        %   params.muDist     (preferred distance in distance units)

        sigmaDist = params.sigmaDist;
        gaussDist = exp(-0.5 * ((dist - params.muDist)^2) / (sigmaDist^2));
        
        val = val+gaussDist;

    case 'angle_vonmises_distance_gaussian'
        % A 2D "von Mises × Gaussian" in (angle × distance) space:
        %   - Von Mises for the angle difference
        %   - Gaussian for the radial distance
        %
        % Required params:
        %   params.muAngle    (preferred angle in degrees)
        %   params.kappaAngle (Von Mises 'concentration' parameter: the larger it is, the more peaked the distribution)
        %   params.muDist     (preferred distance in distance units)
        %   params.sigmaDist  (width in distance units)
    
        % Convert angleDiff (in degrees) to radians
        angleDiffRad = deg2rad(angleDiff);
    
        % Convert muAngle (in degrees) to radians
        muAngleRad = deg2rad(params.muAngle);
    
        % Von Mises concentration parameter
        kappa = params.kappaAngle;
    
        % -- 1) Evaluate the Von Mises for the angle difference
        %
        % Standard Von Mises form (normalized):
        %   VM(θ | μ, κ) = [ exp( κ * cos(θ - μ) ) ] / [2π I₀(κ)]
        % 
        % If you don’t strictly need normalization (e.g. only shape matters),
        % you could omit dividing by (2π I₀(κ)).
    
        vmAngle = exp( kappa * cos(angleDiffRad - muAngleRad) );
    
        % -- 2) Evaluate the Gaussian for distance
        sigmaDist = params.sigmaDist;
        gaussDist = exp( -0.5 * ( (dist - params.muDist)^2 ) ...
                                / (sigmaDist^2) );
    
        % -- 3) Combine them multiplicatively
        val = val + vmAngle * gaussDist;
        
    case 'single_2d_gaussian'
        % A 2D Gaussian centered at an egocentric location.
        %
        % Required params:
        %   params.sigma   (width of the Gaussian in distance units)
        %   params.egocenter   (center of the Gaussian in [x, y] units)

        sigma = params.sigma;
        egoCenter = params.egocenter;  % [cx, cy] in the animal's reference frame

        % Convert heading to radians
        headingRad = deg2rad(head_direction);

        % Construct the rotation matrix for the heading
        R = [ cos(headingRad), -sin(headingRad);
            sin(headingRad),  cos(headingRad) ];

        % Rotate egocenter into global coordinates
        offsetGlobal = R * egoCenter(:);  % ensures it's a column vector

        % Actual center in global space (allocentric)
        centerGlobal = position(:) + offsetGlobal;

        % Now compute distance from centerGlobal to target
        dx_c = targets(i,1) - centerGlobal(1);
        dy_c = targets(i,2) - centerGlobal(2);
        dist_c = sqrt(dx_c^2 + dy_c^2);

        val = val+exp(-0.5 * (dist_c^2) / (sigma^2));
    
    case 'stretched_2d_gaussian'
        % A 2D elliptical Gaussian centered at an egocentric location.
        %
        % Required params (example):
        %   params.sigma_major     (semi-major axis)
        %   params.sigma_minor     (semi-minor axis)
        %   params.ellipse_angle   (orientation of the ellipse, degrees)
        %   params.egocenter       (center of the Gaussian in [x, y] units)

        sigmaX  = params.sigma_major;
        sigmaY  = sigmaX*params.ratio;
        egoCenter = params.egocenter;  % [cx, cy] in the animal's reference frame
        ellipseAngleDeg = mod(head_direction+params.ellipse_angle+180,360)-180;
        % Convert heading to radians (if you need to rotate the center by heading)
        headingRad = deg2rad(head_direction);

        % Construct the rotation matrix for the heading 
        % (used to place the "egocenter" in global coords)
        R_head = [ cos(headingRad), -sin(headingRad);
                sin(headingRad),  cos(headingRad) ];

        % Rotate egocenter into global coordinates
        offsetGlobal = R_head * egoCenter(:);  % ensures it's a column vector

        % Actual center in global space (allocentric)
        centerGlobal = position(:) + offsetGlobal;

        % Now compute the vector from centerGlobal to target (dx_c, dy_c)
        dx_c = targets(i,1) - centerGlobal(1);
        dy_c = targets(i,2) - centerGlobal(2);

        % Convert the ellipse orientation angle to radians
        ellipseAngleRad = deg2rad(ellipseAngleDeg);

        % The rotation matrix that orients the ellipse in the global frame 
        R_ellipse = [  cos(ellipseAngleRad),  sin(ellipseAngleRad);
                    -sin(ellipseAngleRad),  cos(ellipseAngleRad) ];

        % Transform (dx_c, dy_c) into the ellipse's local coordinate system.
        localCoords = R_ellipse * [dx_c; dy_c];
        xLocal = localCoords(1);
        yLocal = localCoords(2);

        % Elliptical exponent:
        %   exp(-0.5 * ( (xLocal^2 / sigmaX^2) + (yLocal^2 / sigmaY^2) ))
        val = val + exp( -0.5 * ( (xLocal^2)/(sigmaX^2) + (yLocal^2)/(sigmaY^2) ) );


    case 'bounded_gaussian'
        % A "fan-shaped" Gaussian: 
        %   - There's a preferred angle (e.g. 0°).
        %   - There's an allowed half-width (e.g. ±30°).
        %   - If target is within that angle range, apply a distance-based Gaussian.
        %   - If target is outside the angle range, it contributes 0.
        %
        % Required params:
        %   params.preferredAngle (preferred angle in degrees)
        %   params.halfAngle (angular half-width)
        %   params.sigmaDist (distance-based falloff)
        %   params.muDist (preferred distance in distance units)

        halfAngle = params.halfWidth;
        sigmaDist = params.sigmaDist;
        angleDeviation = angleDiff - params.prefAngle;
        angleDeviation = mod(angleDeviation + 180, 360) - 180;  % wrap to [-180, 180]

        if abs(angleDeviation) <= halfAngle
            % Within the bounding fan, use a distance-based Gaussian
            val = val+exp(-0.5 * ((dist - params.muDist)^2) / (sigmaDist^2));
        else
            % Outside the fan
            val = val+0;
        end
    
    case 'dog'
        % Difference of Gaussians (DoG) RF:
        %   - Two Gaussian functions with different widths (σ1, σ2).
        %   - The DoG is the difference between these two Gaussians.
        %
        % Required params:
        %   params.ratio (ratio of σ1/σ2)
        %   params.sigma2 (width of the second Gaussian)
        %   params.a2 (amplitude of the second Gaussian relative to the first)
        %   params.egocenter   (center of the Gaussian in [x, y] units)
        egoCenter = params.egocenter;  % [cx, cy] in the animal's reference frame

        % Convert heading to radians
        headingRad = deg2rad(head_direction);

        % Construct the rotation matrix for the heading
        R = [ cos(headingRad), -sin(headingRad);
            sin(headingRad),  cos(headingRad) ];

        % Rotate egocenter into global coordinates
        offsetGlobal = R * egoCenter(:);  % ensures it's a column vector

        % Actual center in global space (allocentric)
        centerGlobal = position(:) + offsetGlobal;

        % Now compute distance from centerGlobal to target
        dx_c = targets(i,1) - centerGlobal(1);
        dy_c = targets(i,2) - centerGlobal(2);
        dist_c = sqrt(dx_c^2 + dy_c^2);
        sigma2 = params.sigma2;
        sigma1 = sigma2 * params.ratio;  % first Gaussian is narrower
        amplitude2 = params.a2;

        gauss1 = exp(-0.5 * ((dist_c^2) / (sigma1^2)));
        gauss2 = amplitude2 * exp(-0.5 * ((dist_c^2) / (sigma2^2)));

        val = val + gauss1 - gauss2;

    otherwise
        error('Unknown rftype: %s', rftype);
end
end
end
