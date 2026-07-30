function spec = getRfSpec(receptive_field_type)
%GETRFSPEC Return shape-parameter names and bounds for one RF type.

type = lower(char(string(receptive_field_type)));

switch type
    case 'angle_distance_gaussian'
        paramNames = {'sigmaAngle', 'muAngle', 'sigmaDist', 'muDist'};
        lb = [5, 0, 1, 0];
        ub = [90, 360, 40, 40];
    case 'angle_gaussian'
        paramNames = {'sigmaAngle', 'muAngle'};
        lb = [5, 0];
        ub = [90, 360];
    case 'distance_gaussian'
        paramNames = {'sigmaDist', 'muDist'};
        lb = [1, 0];
        ub = [40, 40];
    case 'single_2d_gaussian'
        paramNames = {'sigma', 'egocenter_x', 'egocenter_y'};
        lb = [0.1, -30, -30];
        ub = [30, 30, 30];
    case 'stretched_2d_gaussian'
        paramNames = {'sigma_major', 'ratio', 'ellipse_angle', 'egocenter_x', 'egocenter_y'};
        lb = [1, 0.01, 0, -30, -30];
        ub = [60, 1, 180, 30, 30];
    case 'dog'
        paramNames = {'ratio', 'sigma2', 'a2', 'egocenter_x', 'egocenter_y'};
        lb = [0.01, 1, 0, -30, -30];
        ub = [1, 40, 1, 30, 30];
    case 'bounded_gaussian'
        paramNames = {'prefAngle', 'halfWidth', 'sigmaDist', 'muDist'};
        lb = [0, 0, 1, 0];
        ub = [360, 90, 40, 40];
    otherwise
        error('Unknown receptive_field_type: %s', receptive_field_type);
end

spec = struct();
spec.type = type;
spec.paramNames = paramNames;
spec.lb = lb;
spec.ub = ub;
spec.A = [];
spec.b = [];
spec.nParams = numel(paramNames);
end
