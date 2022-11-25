function titleStr = plotCoeffLinearModel(lm,sigColorCode)
% takes a fitted linear model lm and displays the coefficients as a bar
% chart, along with confidence intervals and color coding for p-values. The
% order of the color code is in order of increasing significance.
%% preliminary
if nargin==1
    sigColorCode = signifColorCode(lm.Coefficients.pValue(2:end),...
                        ["lime","yellow","orange","red"]);
end
%% code
CI = [lm.Coefficients.Lower(2:end),lm.Coefficients.Upper(2:end)];

barText = categorical(lm.CoefficientNames(2:end));
barText = reordercats(barText,lm.CoefficientNames(2:end));

bar(barText,lm.Coefficients.Estimate(2:end))
plotCIforEachCat(1:height(CI), CI', sigColorCode)
xticklabels(convert2LatexFormat(lm.CoefficientNames(2:end)))
titleStr = sprintf('n observations = %g',lm.NumObservations);
title(sprintf('n observations = %g',lm.NumObservations))

end
