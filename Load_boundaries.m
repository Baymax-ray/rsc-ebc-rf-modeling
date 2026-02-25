function boundaries_all = Load_boundaries(path)
%LOAD_BOUNDARIES Load boundaries from boundaries.csv in the specified format.
%
% boundaries_all = LOAD_BOUNDARIES(path) reads the file 'boundaries.csv' located at
% the given path. The file should have the columns:
%
% boundary_type appear_timestamp disappear_timestamp x1 y1 x2 y2 cx cy radius start_angle end_angle
%
% - boundary_type: string, 'segment' or 'arc'
% - appear_timestamp, disappear_timestamp: numeric, indicating when boundary appears and disappears，-1 means it is always there
% - For segments: x1, y1, x2, y2 are coordinates of endpoints; arc columns are NaN
% - For arcs: cx, cy, radius, start_angle, end_angle define the arc; segment columns are NaN
%
% Returns:
% boundaries_all: a struct array with fields:
%   .type (string)
%   .appear_timestamp (double) in ms
%   .disappear_timestamp (double) in ms
%   .x1, .y1, .x2, .y2 (double) in cm
%   .cx, .cy, .radius, .start_angle, .end_angle (double) in cm and degree

    % Construct the full path to boundaries.txt
    filePath = fullfile(path, 'boundaries.csv');

    % Set import options (assuming tab-delimited or space-delimited file)
    opts = detectImportOptions(filePath, 'FileType', 'text');
    % Ensure boundary_type is read as text (string)
    opts = setvartype(opts, 'boundary_type', 'string');

    % Read the file into a table
    T = readtable(filePath, opts);

    % Check expected variables (optional)
    expectedVars = ["boundary_type","appear_timestamp","disappear_timestamp", ...
                    "x1","y1","x2","y2","cx","cy","radius","start_angle","end_angle"];
    missingVars = setdiff(expectedVars, T.Properties.VariableNames);
    if ~isempty(missingVars)
        error('Missing expected columns: %s', strjoin(missingVars, ', '));
    end

    % Preallocate struct array
    n = height(T);
    boundaries_all = repmat(struct('type', "", 'appear_timestamp', [], 'disappear_timestamp', [], ...
                                   'x1', [], 'y1', [], 'x2', [], 'y2', [], ...
                                   'cx', [], 'cy', [], 'radius', [], 'start_angle', [], 'end_angle', []), n, 1);

    % Populate the struct array
    for i = 1:n
        boundaries_all(i).type = T.boundary_type(i);
        boundaries_all(i).appear_timestamp = T.appear_timestamp(i);
        boundaries_all(i).disappear_timestamp = T.disappear_timestamp(i);
        boundaries_all(i).x1 = T.x1(i);
        boundaries_all(i).y1 = T.y1(i);
        boundaries_all(i).x2 = T.x2(i);
        boundaries_all(i).y2 = T.y2(i);
        boundaries_all(i).cx = T.cx(i);
        boundaries_all(i).cy = T.cy(i);
        boundaries_all(i).radius = T.radius(i);
        boundaries_all(i).start_angle = T.start_angle(i);
        boundaries_all(i).end_angle = T.end_angle(i);
    end
end
