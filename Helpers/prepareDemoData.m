function demoData = prepareDemoData(config)
%PREPAREDEMODATA Load, validate, and ray-trace the example Patrick neuron.

trackingDir = fullfile(config.dataDir, 'tracking_data');
spikeDir = fullfile(config.dataDir, 'spike_timestamps');
trackingFiles = dir(fullfile(trackingDir, '*.txt'));
spikeFiles = dir(fullfile(spikeDir, '*.txt'));
if numel(trackingFiles) ~= 1
    error('Expected exactly one tracking file in %s; found %d.', ...
        trackingDir, numel(trackingFiles));
end
if numel(spikeFiles) ~= 1
    error('Expected exactly one spike file in %s; found %d.', ...
        spikeDir, numel(spikeFiles));
end

trackingFile = fullfile(trackingFiles(1).folder, trackingFiles(1).name);
[~, neuronName] = fileparts(trackingFiles(1).name);
[~, spikeName] = fileparts(spikeFiles(1).name);
if ~strcmp(neuronName, spikeName)
    error('Tracking and spike files must have the same base name.');
end
spikeFile = fullfile(spikeDir, [neuronName '.txt']);
if ~isfile(spikeFile)
    error('Matching spike file not found: %s', spikeFile);
end

boundaryFile = fullfile(config.dataDir, 'boundaries.CSV');
if ~isfile(boundaryFile)
    boundaryFile = fullfile(config.dataDir, 'boundaries.csv');
end
if ~isfile(boundaryFile)
    error('Boundary file not found in %s.', config.dataDir);
end

tracking = readmatrix(trackingFile);
if size(tracking, 2) ~= 4 || isempty(tracking)
    error('Tracking data must be a nonempty four-column numeric matrix.');
end
if any(~isfinite(tracking(:, 1))) || any(diff(tracking(:, 1)) <= 0)
    error('Tracking timestamps must be finite and strictly increasing.');
end

spikeTimes = readmatrix(spikeFile);
spikeTimes = spikeTimes(:);
if isempty(spikeTimes) || any(~isfinite(spikeTimes))
    error('Spike timestamps must be a nonempty finite numeric vector.');
end

timestamp = tracking(:, 1);
xPosition = tracking(:, 2);
yPosition = tracking(:, 3);
headDirection = tracking(:, 4);

nearestFrame = interp1(timestamp, (1:numel(timestamp))', ...
    spikeTimes, 'nearest', 'extrap');
fired = false(size(timestamp));
firedCount = accumarray(nearestFrame, 1, [numel(timestamp), 1]);
fired(nearestFrame) = true;

invalid = ~isfinite(xPosition) | ~isfinite(yPosition) | ...
    ~isfinite(headDirection) | xPosition < 0 | xPosition > 120 | ...
    yPosition < 0 | yPosition > 120;
timestamp(invalid) = [];
xPosition(invalid) = [];
yPosition(invalid) = [];
headDirection(invalid) = [];
fired(invalid) = [];
firedCount(invalid) = [];
headDirection = mod(headDirection, 360);
fired = fired(:) ~= 0;

boundaries = loadBoundaries(boundaryFile);
sparseness = mean(fired);
fprintf('Prepared %s: %d frames, %d spike frames, sparseness %.6f\n', ...
    neuronName, numel(fired), sum(fired), sparseness);

[egoClosest, ~, ~] = traceClosestBoundaryPoints( ...
    xPosition, yPosition, headDirection, boundaries, ...
    config.numClosest, config.rayAngles, config.useParallel);
folds = makeTimeBlockedFolds(numel(fired), config.nFolds);
[foldNaiveLogLikelihood, naiveCvLogLikelihood] = ...
    computeFoldNaiveLogLikelihood(fired, folds);

demoData = struct();
demoData.name = neuronName;
demoData.timestamp = timestamp;
demoData.xPosition = xPosition;
demoData.yPosition = yPosition;
demoData.headDirection = headDirection;
demoData.fired = fired;
demoData.firedCount = firedCount;
demoData.sparseness = sparseness;
demoData.frame = numel(fired);
demoData.boundaries = boundaries;
demoData.rayAngles = config.rayAngles;
demoData.numClosest = config.numClosest;
demoData.egoClosest = egoClosest;
demoData.folds = folds;
demoData.foldNaiveLogLikelihood = foldNaiveLogLikelihood;
demoData.naiveCvLogLikelihood = naiveCvLogLikelihood;
demoData.fps = config.fps;
demoData.trackingFile = trackingFiles(1).name;
demoData.spikeFile = [neuronName '.txt'];
demoData.boundaryFile = char(string(java.io.File(boundaryFile).getName()));
end

function boundaries = loadBoundaries(boundaryFile)
T = readtable(boundaryFile, 'TextType', 'string');
expected = ["boundary_type", "appear_timestamp", "disappear_timestamp", ...
    "x1", "y1", "x2", "y2", "cx", "cy", "radius", ...
    "start_angle", "end_angle"];
missing = setdiff(expected, string(T.Properties.VariableNames));
if ~isempty(missing)
    error('Boundary file is missing columns: %s', strjoin(missing, ', '));
end
if any(T.boundary_type ~= "segment")
    error('The self-contained Demo currently supports segment boundaries only.');
end

template = struct('type', "", 'appear_timestamp', [], ...
    'disappear_timestamp', [], 'x1', [], 'y1', [], 'x2', [], 'y2', [], ...
    'cx', [], 'cy', [], 'radius', [], 'start_angle', [], 'end_angle', []);
boundaries = repmat(template, height(T), 1);
for i = 1:height(T)
    boundaries(i).type = T.boundary_type(i);
    boundaries(i).appear_timestamp = T.appear_timestamp(i);
    boundaries(i).disappear_timestamp = T.disappear_timestamp(i);
    boundaries(i).x1 = T.x1(i);
    boundaries(i).y1 = T.y1(i);
    boundaries(i).x2 = T.x2(i);
    boundaries(i).y2 = T.y2(i);
    boundaries(i).cx = T.cx(i);
    boundaries(i).cy = T.cy(i);
    boundaries(i).radius = T.radius(i);
    boundaries(i).start_angle = T.start_angle(i);
    boundaries(i).end_angle = T.end_angle(i);
end
end
