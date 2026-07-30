function [colorLimits, colorMap, isConstant] = ...
    getSignedRfColorSettings(RF)
%GETSIGNEDRFCOLORSETTINGS Build an actual-range blue-white-red color map.

lo = min(RF(:));
hi = max(RF(:));
isConstant = lo == hi;

if isConstant
    delta = max(eps(max(1, abs(lo))), 1e-15);
    colorLimits = [lo - delta, hi + delta];
    if lo == 0
        colorMap = ones(257, 3);
    else
        colorMap = signedColorMap(colorLimits(1), colorLimits(2), 257);
    end
else
    colorLimits = [lo, hi];
    colorMap = signedColorMap(lo, hi, 257);
end
end

function colorMap = signedColorMap(lo, hi, nColors)
blue = [33, 102, 172] / 255;
white = [1, 1, 1];
red = [178, 24, 43] / 255;
values = linspace(lo, hi, nColors).';
colorMap = zeros(nColors, 3);

if lo < 0 && hi > 0
    negative = values <= 0;
    positive = ~negative;
    negativeScale = (values(negative) - lo) / (0 - lo);
    positiveScale = values(positive) / hi;
    colorMap(negative, :) = blue + negativeScale .* (white - blue);
    colorMap(positive, :) = white + positiveScale .* (red - white);
    [~, zeroIdx] = min(abs(values));
    colorMap(zeroIdx, :) = white;
elseif hi <= 0
    scale = (values - lo) / (hi - lo);
    colorMap = blue + scale .* (white - blue);
else
    scale = (values - lo) / (hi - lo);
    colorMap = white + scale .* (red - white);
end

colorMap = min(max(colorMap, 0), 1);
end
