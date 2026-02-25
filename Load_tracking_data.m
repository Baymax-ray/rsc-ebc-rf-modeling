function [timestamp, x_position, y_position, head_direction] = Load_tracking_data(path,name)
%LOAD_TRACKING_DATA Load tracking data from a text file.
%
% [timestamp, x_position, y_position, head_direction] = LOAD_TRACKING_DATA(path)
% loads the tracking data located in `path/tracking_data.txt`. 
%
% The file should have four columns:
%   1. timestamp (numeric)
%   2. x_position (numeric)
%   3. y_position (numeric)
%   4. head_direction (numeric, e.g., angle in degrees)
%
% Inputs:
%   path (char or string): Path to the directory containing tracking_data.txt
%
% Outputs:
%   timestamp (double vector) in ms
%   x_position (double vector) in cm
%   y_position (double vector) in cm
%   head_direction (double vector) in degree

    % Construct full file path to tracking_data.txt
    filePath = fullfile(path, name);
    
    % Read the data
    data = readmatrix(filePath);
    
    % Extract columns
    timestamp = data(:, 1);
    x_position = data(:, 2);
    y_position = data(:, 3);
    head_direction = data(:, 4);
end
