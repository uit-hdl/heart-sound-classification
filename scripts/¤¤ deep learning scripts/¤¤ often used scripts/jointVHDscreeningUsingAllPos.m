% In this script I do the following:

% screen for VHD using predicted murmur grade of all positions, and 
% estimate performance metrics using cross-validation. I also investigate
% missed cases. For performance on each position, see estSnSpBalancedWay.

% explore difference in symptom prevalence in detected and missed cases
% (includes hypothesis test)

% investigate ability to detect significant cases (symptomatic
% regurgitation or stenosis), both for each disease separately, or
% screening for all cases jointly. Performance metrics are obtained.

%% estimate PREDICTIVE power of MAX-MURMUR-PREDICTION algorithm
% ¤¤ SET NUMBER OF CV-SPLITS ¤¤
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

%% define container variables
NrowsTot = height(HSdata);
predMatrix   = zeros(NrowsTot,4);
targetMatrix = zeros(NrowsTot,4);
sympAndSickMat = zeros(NrowsTot,4);

%% Estimate performance metrics using cross validation results
Nsplits = 8;
VHDnames = {'AR','MR','AS','MS'};
NrowsTot = height(HSdata);

% define storage variables:
sn = zeros(Nsplits,1);
sp = zeros(Nsplits,1);
ac = zeros(Nsplits,1);
AUCmat = zeros(Nsplits,1);
predictions = cell(Nsplits,1);
activations = cell(Nsplits,1);
Ytrue = cell(Nsplits,1);
JvalTot = cell2mat(CVresults.valTot.J);
HSdata.costumTarget = and(HSdata.MSgrade>0,HSdata.anginaOrDyspnea==0);

% allVHDpredMatrix = zeros(height(HSdata),4);
% plot settings:
plotTrain = false;
plotVal   = false;
plotValTot = true;

% close all
for i=1:Nsplits
    
%     targetType = 'sigSymptVHD';
%     targetType = 'sigAsymptVHD';
%     targetType = 'ARsigSympt'
%     targetType = 'sigVHD31';
    targetType = 'AS';
%     targetType = 'costumTarget';
    predWithMaxGrade = true*(targetType~="AS");
%     predWithMaxGrade = true;
    if targetType(2)=='R'
        classThr = 3;
    else
        classThr = 1;
    end

    I_disease = disease2index(targetType);
    if length(targetType)>2 %#ok<*BDSCI>
        targetVar = targetType;
    else
        targetVar = sprintf('%sgrade',targetType);
    end
    
    ActMatVal = getZeroPaddedActivMatrix(CVresults.val.activations(i,:),...
                                         CVresults.val.J(i,:),NrowsTot);
    ActMatTrain = getZeroPaddedActivMatrix(CVresults.train.activations(i,:),...
                                           CVresults.train.J(i,:),NrowsTot);
    
    plotVal = false;
    ActMat = ActMatTrain + ActMatVal;
    minSn = 0.5;
    minSp = 0.0;
    [activ,u0,Ytarget,Ypred,AUC,glm] = get_sigVHDactivations(ActMat,HSdata,...
                                            CVresults.trainTot.I{i},...
                                            CVresults.valTot.I{i},...
                                            targetType,classThr,...
                                            plotVal,minSn,minSp,predWithMaxGrade);
    % *** store results ***                    
    predictions{i} = Ypred.val;
    activations{i} = activ.val;
    Ytrue{i}       = Ytarget.val;
    
    AUCmat(i) = AUC;
    ac(i) = mean((activ.val>=u0)==Ytarget.val);
    sn(i) = condProb(Ypred.val,Ytarget.val);
    sp(i) = condProb(~Ypred.val,~Ytarget.val);
    
    targetMatrix(CVresults.valTot.I{i},disease2index(targetType)) = Ytarget.val;
    predMatrix(CVresults.valTot.I{i},disease2index(targetType))   = Ypred.val;
    
end

predictions = cell2mat(predictions);
predictionsPadded = predMatrix(:,I_disease);
activations = cell2mat(activations);
YtargetTot = HSdata.(targetVar)(JvalTot)>=classThr;
% YtargetTot = HSdata.ARsigSympt(JvalTot)>=1;
% getAUC(HSdata.ARsigSympt(JvalTot)>=1,activations)
% getAUC(HSdata.MRsigSympt(JvalTot)>=1,activations)

% get table with statistics on SN, SP, and AUC:
T1 = getPredPerf_aucSnSp(predictions,YtargetTot,AUCmat,1)

% *** investigate missed cases ***
S.(targetType) = succAndFailureAnalysis(predictions,YtargetTot,JvalTot,HSdata,...
                {'maxMeanMurGrade','avmeanpg','avarea','anginaOrDyspnea',...
                'sigVHD31','ASgrade','MSgrade','ARgrade','MRgrade',...
                'ARsigSympt','MRsigSympt','sigSymptVHD'});
YtargetTot = HSdata.(targetVar)(JvalTot)>=classThr;

if plotValTot
    hold on
    [~,~,~,~,~,~] = performanceSummaryNeurNet([],YtargetTot,activations,...
                              [],[],plotValTot);
    
    plotSnSpSamples = true;
    if plotSnSpSamples
        colors = {'g','b'};
        thrs = {1, 3.5};
        for i=1:2
            sn = condProb(activations>=thrs{i},YtargetTot);
            sp = condProb(activations<thrs{i},~YtargetTot);
            plot(1-sp,sn,'o','MarkerFaceColor',colors{i})
        end
    end
    
end


% *** estimate AUC and prediction accuracy metrics ***
% get AUC:
AUCci   = computeCImeanEst(AUCmat,"tdist");
aucTot = [mean(AUCmat);AUCci(1);AUCci(2)];
% get SN,SP, and AC:
% store in structure:
thrStr = num2name(classThr);
predPerf.murPred.allPos.(targetType).(thrStr).T = T1;
predPerf.murPred.allPos.(targetType).(thrStr).auc = AUCmat;
predPerf.murPred.allPos.(targetType).(thrStr).T %#ok<*NOPTS>

if targetType(1:2)=='si'
    S.sigVHD31.values.ARsigSympt
end
%%
S.MR.values.MRsigSympt
%% difference in symptomaticness; detected vs missed cases
% Are symptoms more common amongst the true positive AR predictions than
% the false
JvalTot = unionIterated(CVresults.valTot.J);

target = 'AR';
symptomVec = (HSdata.dyspneaCalmlyFlat + HSdata.dyspneaRest + HSdata.angina)>0;
sym.(target).pop = computeCImeanEst(symptomVec,"2")*100;
sym.(target).missedPos = computeCImeanEst(symptomVec(S.(target).J0.falseNeg),"2")*100 %#ok<*NOPTS>
sym.(target).caughtPos = computeCImeanEst(symptomVec(S.(target).J0.truePos),"2")*100
sym.(target).missedNeg = computeCImeanEst(symptomVec(S.(target).J0.falseNeg),"2")*100
sym.(target).allTruePos = computeCImeanEst(...
                symptomVec(HSdata.(sprintf('%sgrade',target))>=classThr),"2")*100;
            
Ivhd = disease2index(target);
IconfirmedVHDandSymp = sympAndSickMat(:,Ivhd);
Ypred = predMatrix(:,Ivhd);
condProb(Ypred(JvalTot),IconfirmedVHDandSymp(JvalTot))
condProb(~Ypred(JvalTot),~IconfirmedVHDandSymp(JvalTot))

sym.(target)

close all
plotCIforEachCat([1,2,3],[sym.(target).pop',...
                          sym.(target).missedPos',...
                          sym.(target).caughtPos'])
%% hypotesis test: are positive predictions more symptomatic?
target = 'AR'
% p1 = symptom prevalence amongst false negative predictions
% p2 = symptom prevalence amongst true positive predictions
p1 = sym.(target).missedPos(2)/100;
p2 = sym.(target).caughtPos(2)/100;
n1 = numel(S.(target).J0.falseNeg);
n2 = numel(S.(target).J0.truePos);
% hypothesis: the proportion of missed AR that has symptoms is higher than
% the proportion of the correctly identified cases of AR that has
% symptoms. In other words, the cases that have been caught have a higher
% prevalence of symptoms than the cases that where missed; p2-p1>0
SE = sqrt(p2*(1-p2)/n1 + p1*(1-p1)/n1);
confInt = p2-p1 + [-1 1]*norminv(0.975)*sqrt(p2*(1-p2)/n1 + p1*(1-p1)/n1);
tvalue = p2-p1; 
presentFraction([p1,p2])
pval = min(1-normcdf(tvalue/SE),normcdf(tvalue/SE))*2
%% (individual) prediction of SYMPTOMATIC VHD:
clear S T
VHDsympt = {'angina','dyspneaCalmlyFlat','dyspneaRest'};
sympSigVHD = getIndexSymptomaticVHD(HSdata,[3,3,1,1],VHDsympt);
for i=1:4   
    S.(VHDnames{i}) = succAndFailureAnalysis(predMatrix(JvalTot,i),sympSigVHD(JvalTot,i),JvalTot,HSdata);
end

myTables.indiv.symp = table((S.AR.summary),...
                            (S.MR.summary),...
                            (S.AS.summary),...
                            (S.MS.summary),'v',VHDnames);

myTables.indiv.symp
%% (joint) screening of SYMPTOMATIC VHD:
clear T S
VHDsympt = {'angina','dyspneaCalmlyFlat','dyspneaRest'};
sympSigVHD = getIndexSymptomaticVHD(HSdata,[3,3,1,1],VHDsympt);

% joint predictions
predAtleastOne = sum(predMatrix,2)>0; %#ok<*NASGU>
YsigSickAtleastOne = sum(sympAndSickMat,2)>0;

Ypred   = predAtleastOne(JvalTot);
Ytarget = YsigSickAtleastOne(JvalTot);
T.sympSigVHD = getTableOverView_signCases(sympSigVHD,predMatrix,JvalTot)
T.sigGradeVHD = getTableOverView_signCases(targetMatrix,predMatrix,JvalTot)
T.sympSigVHD
%%
% writetable(T.sympSigVHD,'symptSigVHDsn50_grade4RegurgOnly.xlsx')
% writetable(T.sigGradeVHD,'sigVHDsn50_grade4RegurgOnly.xlsx')

% sensitivity and specificity for detecting symSig-VHD
S.VHDscreening.sympsig.any.sn = condProb(Ypred,Ytarget);
S.VHDscreening.sympsig.any.sp = condProb(~Ypred,~Ytarget);
S.VHDscreening.sympsig.any.Ncaught = sum(and(Ypred,Ytarget));
S.VHDscreening.sympsig.any.Nmissed = sum(and(~Ypred,Ytarget));

computeCI = true;
% SN and SP for specific diseases for joint-VHD screening:
[S.VHDscreening.sympsig.AR.sn,ciAR] = condProb(Ypred,sympSigVHD(JvalTot,1),computeCI);
[S.VHDscreening.sympsig.MR.sn,ciMR] = condProb(Ypred,sympSigVHD(JvalTot,2),computeCI);
[S.VHDscreening.sympsig.AS.sn,ciAS] = condProb(Ypred,sympSigVHD(JvalTot,3),computeCI);
[S.VHDscreening.sympsig.MS.sn,ciMS] = condProb(Ypred,sympSigVHD(JvalTot,4),computeCI);

[S.VHDscreening.sympsig.AR.sp,ciAR] = condProb(~Ypred,~sympSigVHD(JvalTot,1),computeCI);
[S.VHDscreening.sympsig.MR.sp,ciMR] = condProb(~Ypred,~sympSigVHD(JvalTot,2),computeCI);
[S.VHDscreening.sympsig.AS.sp,ciAS] = condProb(~Ypred,~sympSigVHD(JvalTot,3),computeCI);
[S.VHDscreening.sympsig.MS.sp,ciMS] = condProb(~Ypred,~sympSigVHD(JvalTot,4),computeCI);

% *** individual predictions ***
% Number of cases detected:
S.VHDscreening.sympsig.AR.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,1)));
S.VHDscreening.sympsig.MR.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,2)));
S.VHDscreening.sympsig.AS.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,3)));
S.VHDscreening.sympsig.MS.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,4)));

% Number of cases missed:
S.VHDscreening.sympsig.AR.Nmissed = sum(and(~Ypred,sympSigVHD(JvalTot,1)));
S.VHDscreening.sympsig.MR.Nmissed = sum(and(~Ypred,sympSigVHD(JvalTot,2)));
S.VHDscreening.sympsig.AS.Nmissed = sum(and(~Ypred,sympSigVHD(JvalTot,3)));
S.VHDscreening.sympsig.MS.Nmissed = sum(and(~Ypred,sympSigVHD(JvalTot,4)));
S.VHDscreening.sympsig.MR

T.sympSig.minSn30 = [struct2table(S.VHDscreening.sympsig.AR);...
    struct2table(S.VHDscreening.sympsig.MR);...
    struct2table(S.VHDscreening.sympsig.AS);...
    struct2table(S.VHDscreening.sympsig.MS);...
    struct2table(S.VHDscreening.sympsig.any)]
T.sympSig.minSn30.Properties.RowNames = {'AR','MR','AS','MS','allVHD'};
T.sympSig.minSn30 = rows2vars(T.sympSig.minSn30)
T.sympSig.minSn30

%% (joint) screening of significant VHD:
YpredSigVHD = sum(predMatrix,2)>0;
YtargetSigVHD = sum(targetMatrix,2)>0;

Ssig = succAndFailureAnalysis(Ypred,Ytarget,JvalTot,HSdata,...
        {'maxMeanMurGrade','avmeanpg','avarea'})

S.jointScreen.gradesig.all.sn = condProb(YpredSigVHD(JvalTot),YtargetSigVHD(JvalTot));
S.jointScreen.gradesig.all.sp = condProb(~YpredSigVHD(JvalTot),~YtargetSigVHD(JvalTot));
S.jointScreen.gradesig.all

S.jointScreen.gradesig.AR.sn = condProb(YpredSigVHD(JvalTot),targetMatrix(JvalTot,1));
S.jointScreen.gradesig.MR.sn = condProb(YpredSigVHD(JvalTot),targetMatrix(JvalTot,2));
S.jointScreen.gradesig.AS.sn = condProb(YpredSigVHD(JvalTot),targetMatrix(JvalTot,3));
S.jointScreen.gradesig.MS.sn = condProb(YpredSigVHD(JvalTot),targetMatrix(JvalTot,4));

S.VHDscreening.sympsig.AR.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,1)));
S.VHDscreening.sympsig.AR.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,1)));
S.VHDscreening.sympsig.AR.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,1)));
S.VHDscreening.sympsig.AR.Ncaught = sum(and(Ypred,sympSigVHD(JvalTot,1)));


                      
