function [activ,u0,Ytarget,Ypred,AUC,glm] = get_ASactivations(ActMat,dataFrame,...
                                                Itrain,Ival,classThr,plotVal,...
                                                minSn,minSp)
% computes activations for prediction of AS. ActMat is a 4-by-Ndata (2124
% for HSdata) zero-padded matrix with murmur-algorithm activations.
% Estimates model parameters of multi-position model using AVMEANPG on the
% training set, and then returns activations for the observations
% corresponding to the index vector Ival. Optional output parameters is the
% optimal decision threshold u0 and Ytarget which contains ground truth
% corresponding to predictions, allowing for easy computation of
% sensitivity and specificity.
if nargin==4
    classThr = 1;
    plotVal = false;
    minSn = 0;
    minSp = 0;
elseif nargin==5
    plotVal = false;
    minSn = 0;
    minSp = 0;
elseif nargin==6
    minSn = 0;
    minSp = 0;
end
X = array2table(ActMat,'v',{'A','P','T','M'});
Yavmpg = array2table(dataFrame.avmeanpg,'v',{'avmpg'});
Yas = array2table(dataFrame.ASgrade>=classThr,'v',{'as'});
noise = array2table(ActMat==0,'v',{'noiseA','noiseP','noiseT','noiseM'});

data = [X,noise,Yavmpg,Yas];
I_avmpgTrain = and(Itrain, ~isnan(data.avmpg));
% *** fit model using training data ***
formula = 'avmpg ~ A:A + P:P + T + (A:A):noiseP + noiseP:noiseA:noiseT:M';
% formula = 'avmpg ~ A:A + P:P + T + P:P:noiseA + (A:A):noiseP + noiseP:noiseA';
glm = fitglm(data(I_avmpgTrain,:),formula);
                    

% *** get activation threshold ***
Ytarget.train = dataFrame.ASgrade(Itrain)>=classThr;
activ.train = glm.feval( data(Itrain,:) );
[~,X,Y,T,~] = performanceSummaryNeurNet([],Ytarget.train,activ.train,...
                                        [],[],false);                          
u0 = getOptimalThr(X.whole,Y.whole,T.whole,minSn,minSp);

% *** predict ***
Ytarget.val = dataFrame.ASgrade(Ival)>=classThr;
activ.val   = glm.feval( data(Ival,:) );
Ypred.val   = activ.val>=u0;

if plotVal
    figure
end
[X,Y,~,AUC] = perfcurve(Ytarget.val, activ.val, true);

if plotVal
    plot(X,Y)
    title(sprintf('AUC=%g',round(AUC*100,2)))
end

end