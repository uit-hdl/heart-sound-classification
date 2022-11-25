% This script is basically the same as 'jointVHDscreeningUsingAllPos', but
% it uses the risk factor model for prediction instead.

%% preliminary 1:
net = "regOver";
if net=="reg"
    load CVresults_netMurRegAllPos_valStop
elseif net=="regOver"
    load CVresults_netMurRegAllPos_valStop_overTrain
elseif net=="G2"
    load CVresults_netMurG2AllPos_valStop
elseif net=="G2Over"
    load CVresults_netMurG2AllPos_valStop_overTrain
end
%% Preliminary 2:
NrowsTot = height(HSdata);
predMat   = zeros(NrowsTot,4);
targetMat = zeros(NrowsTot,4);
sympAndSickMat = zeros(NrowsTot,4);
noSympAndSickMat = zeros(NrowsTot,4);
%%
%#ok<*IFBDUP,*BDSCI>
clear perf
%#ok<*NOPTS>
Nsplits = 8;
NrowsTot = height(HSdata);

% ¤¤ SELECT MURMUR VARIABLE ¤¤
% murVar = 'maxMeanMurGrade';
murVar = 'predMaxMurGrade';
% murVar = 'murGrade2';
% murVar = 'predMurGrade';
HSdata.Xmur = zeros(NrowsTot,1);

% ¤¤ SELECT RISK FACTOR MODEL FOR EACH VHD ¤¤
names.all.var = {murVar,'AGE_T7','PULSESPIRO_T72','SEX_T7',...
                'DYSPNEA_FAST_UPHILL_T7','CHEST_PAIN_FAST_T7','HIGH_BLOOD_PRESSURE_T7',...
                'DIABETES_T7','SMOKE_DAILY_Q2_T7'};
names.all.categorical = {'sex','dyspneaFastUpphill','CHEST_PAIN_FAST_T7',...
                'HIGH_BLOOD_PRESSURE_T7','diabetes','smoke','dyspneaCalmlyFlat'};
names.AR.var = {'Xmur','age','sex','dyspneaFastUpphill','pulseSpiro'};
names.MR.var = {'Xmur','age','pulseSpiro'};
names.AS.var = {'Xmur','sex','sex:Xmur'};
names.MS.var = {'Xmur','age','pulseSpiro'};

% ¤¤ SELECT MINIMUM SN AND SP ¤¤
minSn = .5;
minSp = 0;

% ¤¤ SELECT TARGET ¤¤
targetType = 'AR';
% HSdata.ASPGgrade = HSdata.AVMEANPG_T72>=25;
% ¤¤ SET CLASS THRESHOLD ¤¤
if targetType(2)=='R'
    classThr = 3;
elseif targetType(2)=='S'
    if targetType=="ASPG" 
        classThr = 1;
    else
        classThr = 1;
    end
end
plotVal = false;
plotTot = false;
createSubplot = false;
I_disease = disease2index(targetType);

Inan = isnan(table2array(subtable(HSdata,intersect(names.(targetType).var,...
                                         HSdata.Properties.VariableNames) )));
Inan = sum(Inan,2)>0;
Icomplete = ~Inan;
Jcomplete = find(Icomplete);

% *** storage variables ***
AUCmat = zeros(Nsplits,1);
sn = zeros(Nsplits,1);
sp = zeros(Nsplits,1);
activations    = cell(Nsplits,1);
murPredictions = cell(Nsplits,1);
predictions    = cell(Nsplits,1);
Jusable        = cell(Nsplits,1);
Ytrue          = cell(Nsplits,1);
if ~createSubplot
    close all
end
for i=1:8
    ActMatVal = getZeroPaddedActivMatrix(CVresults.val.activations(i,:),...
                                   CVresults.val.J(i,:),NrowsTot);
    ActMatTrain = getZeroPaddedActivMatrix(CVresults.train.activations(i,:),...
                                   CVresults.train.J(i,:),NrowsTot);
    % get activation matrix (training and validation rows do not overlap):
    ActMat = ActMatVal + ActMatTrain;
    covarNames = names.(targetType).var;
    % ¤¤ CHOOSE IF PLOT EACH VALIDATION ROC-CURVE ¤¤
    plotVal = false;
    
    [activ,u0,Ytarget,Ypred,AUC,glm] = get_riskFacModelActivations(ActMat,...
                            HSdata,CVresults.trainTot.I{i},...
                            CVresults.valTot.I{i},Icomplete,...
                            targetType,murVar,...
                            covarNames,names.all.categorical,...
                            classThr,plotVal,minSn,minSp);
                        
    % *** store results ***
    % $ performance $
    AUCmat(i) = AUC;
    Ytrue{i} = Ytarget.val;
    activations{i} = activ.val;
    predictions{i} = Ypred.val;
    sn(i) = condProb(Ypred.val,Ytarget.val);
    sp(i) = condProb(~Ypred.val,~Ytarget.val);
    
    % $ store target and predictions in full-form matrices $
    targetMat(and(CVresults.valTot.I{i},Icomplete),I_disease) = Ytarget.val;
    predMat(  and(CVresults.valTot.I{i},Icomplete),I_disease) = Ypred.val;

end

% variables for the set of all CV-predictions:
Jval = cell2mat(CVresults.valTot.J);
activations = cell2mat(activations);
predictions = cell2mat(predictions);
Ytrue       = cell2mat(Ytrue);

% table that summarizes info on SN, SP and AUC:
T = getPredPerf_aucSnSp(predictions,Ytrue,AUCmat,1);
thrStr = num2name(classThr)
predPerf.riskFac.(targetType).(thrStr).T = T;
predPerf.riskFac.(targetType).(thrStr).auc = AUCmat;
predPerf.riskFac.(targetType).(thrStr).T

% S contains info that allows inspection into correct and incorrect
% predictions:
S.(targetType) = succAndFailureAnalysis(predictions,Ytrue,Jval,HSdata,...
                        {'maxMeanMurGrade','avmeanpg','avarea'});
if plotTot
    AUC = performanceSummaryNeurNet([],Ytrue,activations,...
                                 numel(activations),[],plotTot);
    title(sprintf('AUC=%g%%',round(AUC.whole,3)*100))
end
compactLinModelPresentation(glm)
%% prediction of significant-symptomatic VHD:
Jpred = unionIterated(CVresults.valTot.J);

% joint predictions
predAtleastOne = sum(predMat,2)>0;
sigSickAtleastOne = sum(sympAndSickMat,2)>0;
sigSickAndNoSymp = sum(noSympAndSickMat,2)>0;

Ypred = predAtleastOne(Jpred);
Ytarget = sigSickAtleastOne(Jpred);

sympVHD.all.sn = condProb(Ypred,Ytarget);
sympVHD.all.sp = condProb(~Ypred,~Ytarget);

% *** individual predictions ***
% how many detected:
sympVHD.AR.Ncaught = sum(and(Ypred,sympAndSickMat(Jpred,1)));
sympVHD.MR.Ncaught = sum(and(Ypred,sympAndSickMat(Jpred,2)));
sympVHD.AS.Ncaught = sum(and(Ypred,sympAndSickMat(Jpred,3)));
sympVHD.MS.Ncaught = sum(and(Ypred,sympAndSickMat(Jpred,4)));

% how many not detected:
sympVHD.AR.Nmissed = sum(and(~Ypred,sympAndSickMat(Jpred,1)));
sympVHD.MR.Nmissed = sum(and(~Ypred,sympAndSickMat(Jpred,2)));
sympVHD.AS.Nmissed = sum(and(~Ypred,sympAndSickMat(Jpred,3)));
sympVHD.MS.Nmissed = sum(and(~Ypred,sympAndSickMat(Jpred,4)));

% how many not detected:
sympVHD.AR.sn = condProb(Ypred,sympAndSickMat(Jpred,1));
sympVHD.MR.sn = condProb(Ypred,sympAndSickMat(Jpred,2));
sympVHD.AS.sn = condProb(Ypred,sympAndSickMat(Jpred,3));
sympVHD.MS.sn = condProb(Ypred,sympAndSickMat(Jpred,4));


sympVHD.all

%% prediction of significant VHD:
YpredSigVHD = sum(predMat,2)>0;
YtargetSigVHD = sum(targetMat,2)>0;
sigVHD.all.sn = condProb(YpredSigVHD(Jpred),YtargetSigVHD(Jpred));
sigVHD.all.sp = condProb(~YpredSigVHD(Jpred),~YtargetSigVHD(Jpred));
sigVHD.all

sigVHD.AR.sn = condProb(YpredSigVHD(Jpred),targetMat(Jpred,1));
sigVHD.MR.sn = condProb(YpredSigVHD(Jpred),targetMat(Jpred,2));
sigVHD.AS.sn = condProb(YpredSigVHD(Jpred),targetMat(Jpred,3));
sigVHD.MS.sn = condProb(YpredSigVHD(Jpred),targetMat(Jpred,4));

sigVHD.AR
sigVHD.MR
%% difference in symptomaticness; detected vs missed cases
Jpred = unionIterated(CVresults.valTot.J);
sym.(targetVarType).pop = computeCImeanEst(symptomVec,"2")*100;
sym.(targetVarType).missedPos = computeCImeanEst(symptomVec(JfalseNeg),"2")*100 %#ok<*NOPTS>
sym.(targetVarType).caughtPos = computeCImeanEst(symptomVec(JtruePos),"2")*100
sym.(targetVarType).missedNeg = computeCImeanEst(symptomVec(JfalseNeg),"2")*100
sym.(targetVarType).allTruePos = computeCImeanEst(symptomVec(HSdata.(targetVar)>=classThr),"2")*100
Ivhd = disease2index(targetVarType);
IconfirmedVHDandSymp = sympAndSickMat(:,Ivhd);
Ypred = predMat(:,Ivhd);
condProb(Ypred(Jpred),IconfirmedVHDandSymp(Jpred))

condProb(~Ypred(Jpred),~IconfirmedVHDandSymp(Jpred))

sym.(targetVarType)

close all
plotCIforEachCat([1,2,3],[sym.(targetVarType).pop',...
                          sym.(targetVarType).missedPos',...
                          sym.(targetVarType).caughtPos'])
                      
%% hypotesis test
p1 = sym.(targetVarType).missedPos(2)/100;
p2 = sym.(targetVarType).caughtPos(2)/100;
n1 = numel(JfalseNeg);
n2 = numel(JtruePos);

% hypothesis: the proportion of missed AR that has symptoms is higher than
% the proportion of the correctly identified cases of AR that has
% symptoms. In other words, the cases that have been caught have a higher
% prevalence of symptoms than the cases that where missed; p2-p1>0
SE = sqrt(p2*(1-p2)/n1 + p1*(1-p1)/n1);
confInt = p2-p1 + [-1 1]*norminv(0.975)*sqrt(p2*(1-p2)/n1 + p1*(1-p1)/n1);
tvalue = p2-p1;
1-normcdf(tvalue/SE)