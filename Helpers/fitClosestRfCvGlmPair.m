function [freeResult, posResult] = fitClosestRfCvGlmPair(Egorelative_closest, fired, ...
    receptive_field_type, numClosest, folds, foldNaiveLogLikelihood, fitOptions)
%FITCLOSESTRFCVGLMPAIR Fit glm_free first, then reuse valid glm_pos folds.

if nargin < 7 || isempty(fitOptions)
    fitOptions = struct();
end

freeResult = fitClosestRfCvCore(Egorelative_closest, fired, receptive_field_type, ...
    numClosest, 'glm_free', folds, foldNaiveLogLikelihood, fitOptions, []);
posResult = fitClosestRfCvCore(Egorelative_closest, fired, receptive_field_type, ...
    numClosest, 'glm_pos', folds, foldNaiveLogLikelihood, fitOptions, freeResult);
end
