function logLikelihood = logisticLogLikelihood(fired, linearPredictor)
%LOGISTICLOGLIKELIHOOD Bernoulli log likelihood for logistic probabilities.

epsval = 1e-15;
fired = double(fired(:) ~= 0);
linearPredictor = double(linearPredictor(:));

prob = 1 ./ (1 + exp(-linearPredictor));
prob = min(max(prob, epsval), 1 - epsval);

logLikelihood = sum(fired .* log(prob) + (1 - fired) .* log(1 - prob));
end
