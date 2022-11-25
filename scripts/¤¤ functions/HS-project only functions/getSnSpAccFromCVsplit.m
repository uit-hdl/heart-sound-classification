function [perf,Ypred] = getSnSpAccFromCVsplit(Xtrain,Xval,Ytrain,Yval,minSn,minSp)
% Estimates sn, sp, and accuracy on validation set, using decision
% threshold that is estimated using data from vectors Xtrain an Ytrain.
%% preliminary
if nargin==4
    minSn = 0;
    minSp = 0;
elseif nargin==5
    minSp = 0;
end
%% example
% Xtrain = predTrain;
% Xval   = predVal;
% Ytrain = YtrueTrain;
% Yval   = YtrueVal;
% minSensitivity = 0;
%% code:
% *** estimate optimal threshold ***
[X,Y,T] = perfcurve(Ytrain,Xtrain,true);
u0      = getOptimalThr(X,Y,T,minSn,minSp);
Ypred = Xval>=u0;

sn = condProb(Ypred,Yval);
sp = condProb(~Ypred,~Yval);
ac = prob(Ypred==Yval);

perf = [sn,sp,ac];

end