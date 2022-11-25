function [perf,Ypred] = getPredMetricsAndPred(Xtrain,Xval,Ytrain,Yval,minSensitivity)
% Computes
% Xtrain and Xval are vectors containing activations for training and
% validation set respectively, and similarly Ytrain and Yval are a
% ground-truth logical index vectors of their respective sets. Returns
% sensitivity, specificity and accuracy with confidence intervals computed
% using exact methods (assumes only that probability of success is
% independent) based on the PDF of the binomial distribution.
%% example
% Xtrain = predTrain;
% Xval   = predVal;
% Ytrain = YtrueTrain;
% Yval   = YtrueVal;
% minSensitivity = 0;
%% code:
% *** estimate optimal threshold ***
[X,Y,T] = perfcurve(Ytrain,Xtrain,true);
u0      = getOptimalThr(X,Y,T,minSensitivity);
Ypred = Xval>=u0;

sn = condProb(Ypred,Yval);
sp = condProb(~Ypred,~Yval);
ac = prob(Ypred==Yval);

perf = [sn,sp,ac];

end