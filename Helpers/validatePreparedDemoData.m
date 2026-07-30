function validatePreparedDemoData(demoData, config)
%VALIDATEPREPAREDDEMODATA Ensure cached data match the active configuration.

required = {'name', 'fired', 'frame', 'sparseness', 'egoClosest', ...
    'folds', 'foldNaiveLogLikelihood', 'numClosest', 'rayAngles', 'fps'};
for i = 1:numel(required)
    if ~isfield(demoData, required{i})
        error('Cached demo data are missing field %s.', required{i});
    end
end
if demoData.frame ~= numel(demoData.fired)
    error('Cached frame count does not match fired.');
end
if demoData.numClosest ~= config.numClosest
    error('Cached numClosest does not match the active configuration.');
end
if numel(demoData.folds) ~= config.nFolds
    error('Cached fold count does not match the active configuration.');
end
if ~isequal(demoData.rayAngles, config.rayAngles)
    error('Cached ray angles do not match the active configuration.');
end
if demoData.fps ~= config.fps
    error('Cached fps does not match the active configuration.');
end
if size(demoData.egoClosest, 1) ~= demoData.frame || ...
        size(demoData.egoClosest, 2) ~= config.numClosest
    error('Cached closest-boundary array has an incompatible size.');
end
end
