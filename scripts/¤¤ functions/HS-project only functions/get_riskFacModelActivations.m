function [activ,u0,Ytarget,Ypred,AUC,glm] = get_riskFacModelActivations(ActMat,...
                                                dataFrame,Itrain,Ival,...
                                                Icomplete,targetVarType,murVar,...
                                                namesVar,namesCategorical,...
                                                classThr,plotVal,minSn,minSp)
% computes activations for prediction of VHD specified by targetVarType. A
% a risk-factor model with covariates given by varNames (a character
% array). ActMat is a 4-by-Ndata (2124 for HSdata) zero-padded matrix with
% murmur-algorithm activations. Which murmur-variable is used is indicated
% by the string murVar, which can be of the forms "predMaxMurGrade",
% "predMuri" (i denotes position used), and "maxMeanMurGrade". Note that if
% "predMaxMurGrade" set, then AS will be predicted using the multiposition
% model. Returns activations (reduced form) for the observations
% corresponding to the index vector Ival.

% Optional output parameters is the optimal decision threshold u0, Ytarget
% (ground truth of validation set) and Ypred (validation set predicions),
% allowing for easy computation of sensitivity and specificity.

if nargin==4
    classThr = 1;
    plotVal = false;
    minSn = 0;
    minSp = 0;
elseif nargin==5
    plotVal = false;
    minSn = 0;
    minSp = 0;
end

% *** get murmur-grade activations, and add to dataframe as variable  ***
targetVarName = sprintf('%sgrade',targetVarType);

if murVar=="predMaxMurGrade"
    murActiv = get_sigVHDactivations(ActMat,dataFrame,Itrain,...
                                    Ival,targetVarType,classThr,plotVal,...
                                    minSn,minSp); 
    dataFrame.Xmur(Itrain) = murActiv.train;
    dataFrame.Xmur(Ival)   = murActiv.val;
    
elseif murVar(1:7)=="predMur"
    murActiv = ActMat(:,str2double(murVar(end)));
    dataFrame.Xmur(Itrain) = murActiv(Itrain);
    dataFrame.Xmur(Ival)   = murActiv(Ival);
    
elseif murVar=="maxMeanMurGrade"
    dataFrame.Xmur = dataFrame.(murVar);
end
                           
% add AR, MR, AS or MS as a variable in the data-frame:
dataFrame.(targetVarType) = dataFrame.(targetVarName)>=classThr;

% get the formula:
formula = getLinearModelFormula(namesVar,targetVarType);

% *** fit model using training data ***
glm = fitglm(dataFrame(Itrain,:),formula,...
    'categorical',intersect(namesCategorical,namesVar));%,'distr','binomial');
% *** predict on training set and get activation threshold ***
Ytarget.train = dataFrame.(targetVarName)(and(Itrain,Icomplete))>=classThr;
activ.train = glm.feval(dataFrame(and(Itrain,Icomplete),:));
[~,X,Y,T,~] = performanceSummaryNeurNet([],Ytarget.train,activ.train,...
                                        [],[],false);   

u0 = getOptimalThr(X.whole,Y.whole,T.whole,minSn,minSp);

% *** predict on validation set ***
Ytarget.val = dataFrame.(targetVarName)(and(Ival,Icomplete))>=classThr;
activ.val   = glm.feval(dataFrame(and(Ival,Icomplete),:));
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