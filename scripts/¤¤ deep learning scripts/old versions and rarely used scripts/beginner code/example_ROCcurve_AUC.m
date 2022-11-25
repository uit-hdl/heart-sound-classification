% ROC curve for test data
load fisheriris
% predictors:
pred = meas(51:end,1:2);
% define binary response variable:
resp = (1:100)'>50;  % Versicolor = 0, virginica = 1
mdl = fitglm(pred,resp,'Distribution','binomial','Link','logit');
% get probabilities for belonging to the "positive" class:
p_posClass = mdl.Fitted.Probability;
trueLabels = species(51:end,:);
labelPosClass = 'virginica';

[X,Y,T,AUC] = perfcurve(trueLabels, p_posClass,labelPosClass);
[X,Y,T,AUC] = perfcurve(categorical(trueLabels), p_posClass,'virginica')

% ROC curve for murmur classification
pMur = act(2,:);
pMurTot = reshape(pMur,[Nval,NsegDesired]);
pMurTot = median(pMurTot,2);
[Xseg,Yseg,~,AUCseg] = perfcurve(Yval,pMur,'1');
[Xtot,Ytot,~,AUCtot] = perfcurve(YvalTot,pMurTot,'1');
clf
plot(1-Xseg,Yseg)
hold on
plot(1-Xtot,Ytot)
title(sprintf('AUCseg=%.2g, AUCwhole=%.2g',AUC,AUCtot))
legend({'segments','whole PCG-signal'})


