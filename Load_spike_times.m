function spike_times_real = Load_spike_times(path,name)
%LOAD_SPIKE_TIMES Load spike times from a text file.
%
% spike_times_real = LOAD_SPIKE_TIMES(path) loads the spike times located in 
% `path/spike_times.txt`.
%
% The file should have a single column of numeric values representing the spike times.
%
% Input:
%   path (char or string): Path to the directory containing spike_times.txt
%
% Output:
%   spike_times_real (double vector): Vector of spike times in ms

    % Construct full file path to spike_times.txt
    filePath = fullfile(path, name);
    
    % Read the data (assuming a single column)
    spike_times_real = readmatrix(filePath);
end
