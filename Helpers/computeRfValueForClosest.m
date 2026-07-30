function rfValue = computeRfValueForClosest(receptive_field_type, Egorelative_closest, shapeParams)
%COMPUTERFVALUEFORCLOSEST Sum RF response across closest egocentric points.

param = buildRfParams(receptive_field_type, shapeParams);
[nFrames, nClosest, ~] = size(Egorelative_closest);
pts = reshape(Egorelative_closest, [nFrames * nClosest, 2]);

rfPerPoint = computePointRf(receptive_field_type, pts, param);
rfPerPoint = reshape(rfPerPoint, [nFrames, nClosest]);
rfValue = sum(rfPerPoint, 2);
end

function rfValue = computePointRf(receptive_field_type, data, param)
x = data(:, 1);
y = data(:, 2);

switch lower(char(string(receptive_field_type)))
    case 'single_2d_gaussian'
        distSq = (x - param.egocenter(1)).^2 + (y - param.egocenter(2)).^2;
        rfValue = exp(-distSq / (2 * param.sigma^2));

    case 'angle_distance_gaussian'
        dist = sqrt(x.^2 + y.^2);
        theta = atan2d(y, x);
        diffTheta = mod(theta - param.muAngle + 180, 360) - 180;
        gaussDist = exp(-((dist - param.muDist).^2) / (2 * param.sigmaDist^2));
        gaussAngle = exp(-(diffTheta.^2) / (2 * param.sigmaAngle^2));
        rfValue = gaussDist .* gaussAngle;

    case 'angle_gaussian'
        theta = atan2d(y, x);
        diffTheta = mod(theta - param.muAngle + 180, 360) - 180;
        rfValue = exp(-(diffTheta.^2) / (2 * param.sigmaAngle^2));

    case 'distance_gaussian'
        dist = sqrt(x.^2 + y.^2);
        rfValue = exp(-((dist - param.muDist).^2) / (2 * param.sigmaDist^2));

    case 'stretched_2d_gaussian'
        dx = x - param.egocenter(1);
        dy = y - param.egocenter(2);
        sigmaX = param.sigma_major;
        sigmaY = sigmaX * param.ratio;
        xRot = dx * cosd(param.ellipse_angle) + dy * sind(param.ellipse_angle);
        yRot = -dx * sind(param.ellipse_angle) + dy * cosd(param.ellipse_angle);
        rfValue = exp(-((xRot.^2) / (2 * sigmaX^2) + (yRot.^2) / (2 * sigmaY^2)));

    case 'bounded_gaussian'
        theta = atan2d(y, x);
        diffTheta = mod(theta - param.prefAngle + 180, 360) - 180;
        inFan = abs(diffTheta) < param.halfWidth;
        dist = sqrt(x.^2 + y.^2);
        val = exp(-((dist - param.muDist).^2) / (2 * param.sigmaDist^2));
        rfValue = zeros(size(x));
        rfValue(inFan) = val(inFan);

    case 'dog'
        distSq = (x - param.egocenter(1)).^2 + (y - param.egocenter(2)).^2;
        sigma1 = param.sigma2 * param.ratio;
        gauss1 = exp(-distSq / (2 * sigma1^2));
        gauss2 = param.a2 * exp(-distSq / (2 * param.sigma2^2));
        rfValue = gauss1 - gauss2;

    otherwise
        error('Unknown receptive_field_type: %s', receptive_field_type);
end
end
