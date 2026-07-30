function fit = fitLogisticCalibration(rfValue, fired, calibration)
%FITLOGISTICCALIBRATION Fit alpha/offset for a fixed RF value vector.
%
% glm_free: logit(p) = b0 + b1 * rfValue, alpha = b1, offset = -b0.
% glm_pos: same fit unless b1 < 0, then alpha is constrained to 0.

if nargin < 3 || isempty(calibration)
    calibration = 'glm_free';
end

calibration = lower(char(string(calibration)));
rfValue = double(rfValue(:));
fired = double(fired(:) ~= 0);

valid = isfinite(rfValue) & isfinite(fired);
rfValue = rfValue(valid);
fired = fired(valid);

if isempty(fired)
    error('Cannot fit logistic calibration on empty data.');
end

[beta, method] = fitFreeLogistic(rfValue, fired);
boundaryApplied = false;

switch calibration
    case 'glm_free'
        % use the free fit
    case 'glm_pos'
        if beta(2) < 0
            beta = [logitClipped(mean(fired)); 0];
            method = 'intercept_only_positive_boundary';
            boundaryApplied = true;
        end
    otherwise
        error('Unknown calibration: %s', calibration);
end

linearPredictor = beta(1) + beta(2) * rfValue;

fit = struct();
fit.calibration = calibration;
fit.alpha = beta(2);
fit.offset = -beta(1);
fit.beta = beta(:).';
fit.logLikelihood = logisticLogLikelihood(fired, linearPredictor);
fit.n = numel(fired);
fit.method = method;
fit.boundaryApplied = boundaryApplied;
end

function [beta, method] = fitFreeLogistic(rfValue, fired)
if numel(unique(fired)) < 2 || all(abs(rfValue - rfValue(1)) < 1e-12)
    beta = [logitClipped(mean(fired)); 0];
    method = 'intercept_only';
    return
end

try
    oldWarn = warning('off', 'all');
    cleanupObj = onCleanup(@() warning(oldWarn));
    beta = glmfit(rfValue, fired, 'binomial', 'link', 'logit');
    clear cleanupObj
    method = 'glmfit';
catch
    if exist('cleanupObj', 'var')
        clear cleanupObj
    end
    beta = fitByFminsearch(rfValue, fired);
    method = 'fminsearch';
end

if numel(beta) ~= 2 || any(~isfinite(beta))
    beta = fitByFminsearch(rfValue, fired);
    method = 'fminsearch';
end
end

function beta = fitByFminsearch(rfValue, fired)
init = [logitClipped(mean(fired)); 0];
obj = @(b) -logisticLogLikelihood(fired, b(1) + b(2) * rfValue);
opts = optimset('Display', 'off', 'MaxIter', 1000, 'MaxFunEvals', 5000);
beta = fminsearch(obj, init, opts);
beta = beta(:);
end

function x = logitClipped(p)
epsval = 1e-6;
p = min(max(p, epsval), 1 - epsval);
x = log(p / (1 - p));
end
