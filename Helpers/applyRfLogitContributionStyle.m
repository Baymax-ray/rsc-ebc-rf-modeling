function applyRfLogitContributionStyle(fig, RF, fullFitOffset)
%APPLYRFLOGITCONTRIBUTIONSTYLE Apply the shared RF contribution styling.
% fullFitOffset corresponds to beta in the manuscript.

if any(~isfinite(RF(:)))
    error('RF contains non-finite values.');
end
if ~isnumeric(fullFitOffset) || ~isscalar(fullFitOffset) || ...
        ~isfinite(fullFitOffset)
    error('fullFitOffset must be a finite numeric scalar.');
end

ax = gca;
[colorLimits, colorMap, isConstant] = getSignedRfColorSettings(RF);
colormap(fig, colorMap);
caxis(ax, colorLimits);

cb = colorbar(ax);
cb.Label.String = 'Contribution to firing log-odds (\alpha f_\Theta({\bf b}))';
cb.Label.Interpreter = 'tex';
if isConstant
    cb.Ticks = RF(1);
    cb.TickLabels = {sprintf('%.4g', RF(1))};
end

title(ax, 'Logistic Input Contribution (Egocentric)');
baselineProbability = 1 / (1 + exp(fullFitOffset));
subtitle(ax, sprintf( ...
    'Zero-input baseline: p_0 = %.3g (offset = %.3g)', ...
    baselineProbability, fullFitOffset), 'Interpreter', 'tex');
end
