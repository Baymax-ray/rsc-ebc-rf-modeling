function param = buildRfParams(receptive_field_type, shapeParams)
%BUILDRFPARAMS Convert an RF shape-parameter vector into a struct.

type = lower(char(string(receptive_field_type)));
p = shapeParams(:).';

switch type
    case 'angle_distance_gaussian'
        param.sigmaAngle = p(1);
        param.muAngle = p(2);
        param.sigmaDist = p(3);
        param.muDist = p(4);
    case 'angle_gaussian'
        param.sigmaAngle = p(1);
        param.muAngle = p(2);
    case 'distance_gaussian'
        param.sigmaDist = p(1);
        param.muDist = p(2);
    case 'single_2d_gaussian'
        param.sigma = p(1);
        param.egocenter = [p(2), p(3)];
    case 'stretched_2d_gaussian'
        param.sigma_major = p(1);
        param.ratio = p(2);
        param.ellipse_angle = p(3);
        param.egocenter = [p(4), p(5)];
    case 'dog'
        param.ratio = p(1);
        param.sigma2 = p(2);
        param.a2 = p(3);
        param.egocenter = [p(4), p(5)];
    case 'bounded_gaussian'
        param.prefAngle = p(1);
        param.halfWidth = p(2);
        param.sigmaDist = p(3);
        param.muDist = p(4);
    otherwise
        error('Unknown receptive_field_type: %s', receptive_field_type);
end
end
