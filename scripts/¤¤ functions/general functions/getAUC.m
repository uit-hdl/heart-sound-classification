function [AUC,X,Y,T] = getAUC(Ytrue,Ypred)
% computes the AUC. Assumes Ytrue is logical, and Ypred numerical.
[X,Y,T,AUC] = perfcurve(Ytrue, Ypred, true);
end