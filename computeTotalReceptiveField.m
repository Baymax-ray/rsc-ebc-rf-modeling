function rfValue = computeTotalReceptiveField(receptive_field_type, data, param)
%   Compute the value of a receptive field for each points in data without activation function
%   rfValue = computeReceptiveField(receptive_field_type, data, param)
%
%   INPUTS:
%     receptive_field_type: string specifying RF model, e.g.
%        'single_2d_gaussian', 'angle_distance_gaussian', ...
%     data: Nx2 matrix, each row = (x_i, y_i) location (egocentric coordinates)
%     param: struct containing fields needed by the chosen RF type
%
%   OUTPUT:
%     rfValue: Nx1 vector of RF values at each (x_i, y_i)

    % Separate the x,y coordinates
    x = data(:,1);
    y = data(:,2);

    switch lower(receptive_field_type)

        case 'single_2d_gaussian'
            % Standard 2D Gaussian centered at an egocentric location

            x0     = param.egocenter(1); % center x coordinate
            y0     = param.egocenter(2); % center y coordinate
            sigma  = param.sigma;   % radial std

            % Squared distance from center
            distSq = (x - x0).^2 + (y - y0).^2;
            % Standard circular 2D Gaussian
            rfValue = exp(-distSq/(2*sigma^2));

        case 'angle_distance_gaussian'
            % A 2D Gaussian in "angle × distance" space:
            %   - One Gaussian for the angle difference
            %   - One Gaussian for the radial distance
            % The final RF is the product of these two Gaussians.
            sigmaDist  = param.sigmaDist;
            sigmaAngle = param.sigmaAngle; % (in degrees)
            muDist    = param.muDist;    % (in cm)
            muAngle   = param.muAngle;   % (in degrees)

            % Distance and angle (polar coordinates)
            dist= sqrt(x.^2 + y.^2); % radial distance
            theta = atan2d(y, x);  % angle in degrees
            diffTheta_deg = mod(theta - muAngle + 180, 360) - 180; % angle difference (in degrees)

            gaussDist  = exp(-((dist - muDist).^2)/(2*sigmaDist^2)); % distance Gaussian
            gaussAngle = exp(-((diffTheta_deg).^2)/(2*sigmaAngle^2)); % angle Gaussian

            rfValue = gaussDist .* gaussAngle;

        case 'angle_gaussian'
            % 1D Gaussian in "angle" space:
            %   - One Gaussian for the angle difference
            sigmaAngle = param.sigmaAngle; % (in degrees)
            muAngle   = param.muAngle;   % (in degrees)
            theta = atan2d(y, x);  % angle in degrees
            diffTheta_deg = mod(theta - muAngle + 180, 360) - 180; % angle difference (in degrees)

            rfValue = exp(-((diffTheta_deg).^2)/(2*sigmaAngle^2)); % angle Gaussian

        case 'distance_gaussian'
            % 1D Gaussian in "distance" space:
            %   - One Gaussian for the radial distance
            sigmaDist  = param.sigmaDist; % (in cm)
            muDist    = param.muDist;    % (in cm)

            dist= sqrt(x.^2 + y.^2); % radial distance

            rfValue = exp(-((dist - muDist).^2)/(2*sigmaDist^2)); % distance Gaussian

        case 'stretched_2d_gaussian'
            % 2D elliptical Gaussian centered at an egocentric location
            %   -one major axis (sigmaX) and one minor axis (sigmaY)
            x0     = param.egocenter(1); % center x coordinate
            y0     = param.egocenter(2);
            sigmaX  = param.sigma_major;
            sigmaY  = sigmaX*param.ratio;
            phi    = param.ellipse_angle; % angle of rotation (in degrees)

            % distance to center
            dx = x - x0;
            dy = y - y0;

            % Rotate by -phi to align with principal axes
            x_rot =  dx*cosd(phi) + dy*sind(phi);
            y_rot = -dx*sind(phi) + dy*cosd(phi);

            % Elliptical exponent
            rfValue = exp(-((x_rot.^2)/(2*sigmaX^2) + (y_rot.^2)/(2*sigmaY^2)));

        case 'bounded_gaussian'
            % "Fan-shaped" or "sector-limited" Gaussian, 
            %   - There's a preferred angle (e.g. 0°).
            %   - There's an allowed half-width (e.g. ±30°).
            %   - If target is within that angle range, apply a distance-based Gaussian.
            %   - If target is outside the angle range, it contributes 0.
            muDist    = param.muDist;    % (in cm)
            sigmaDist  = param.sigmaDist; % (in cm)
            prefAngle = param.prefAngle; % (in degrees)
            halfWidth = param.halfWidth; % (in degrees)

            theta= atan2d(y, x);  % angle in degrees
            diffTheta_deg = mod(theta - prefAngle + 180, 360) - 180; % angle difference (in degrees)
            infan = abs(diffTheta_deg) < halfWidth; % inside the fan angle
            
            distSq = x.^2 + y.^2; % squared distance from center
            val=exp(-((sqrt(distSq) - muDist).^2)/(2*sigmaDist^2)); % distance Gaussian
            

            rfValue = zeros(size(x)); % initialize RF value with zeros
            rfValue(infan) = val(infan); % apply distance Gaussian only within the fan angle

        case 'dog'
            % Difference of two Gaussians (DoG) with different sigmas
            %   - Two Gaussian functions with different widths (σ1, σ2).
            %   - The DoG is the difference between these two Gaussians.
            x0     = param.egocenter(1); % center x coordinate
            y0     = param.egocenter(2); % center y coordinate
            sigma2  = param.sigma2;   % radial std of the second Gaussian
            sigma1  = sigma2*param.ratio;   % radial std of the first Gaussian
            a2    = param.a2;    % amplitude of the second Gaussian (the amplitude of the first Gaussian is 1)

            distSq = (x - x0).^2 + (y - y0).^2;
            % First Gaussian (amplitude = 1)
            gauss1 = exp(-distSq/(2*sigma1^2));
            % Second Gaussian (amplitude = a2)
            gauss2 = a2 * exp(-distSq/(2*sigma2^2));

            rfValue = (gauss1 - gauss2);

        otherwise
            error('Unknown receptive_field_type: %s', receptive_field_type);
    end
end
