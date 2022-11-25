% Primary use is to investigate predictive performance of algorithm in each
% position. I mostly only need the first 2 sections, where I gather CV
% performance metrics for prediction of various targets for each position.

%#ok<*NOPTS>
%#ok<*NASGU>
%#ok<*STCMP,*BDSCI,*BDSCA>
%% load data
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

%% Prediction performance for each INDIVIDUAL location
clear AUCmat
% *** estimate sensitivity and specificity for each location ***

% ¤¤ SET NUMBER OF CV-SPLITS ¤¤
Nsplits = 8;
% storage variables:
sn = zeros(Nsplits,4);
sp = zeros(Nsplits,4);
ac = zeros(Nsplits,4); 
AUCmat.val = zeros(Nsplits,5);
AUCmat.train = zeros(Nsplits,4);
predictions  = cell(Nsplits,4);
activations  = cell(Nsplits,4); 
Ytrue        = cell(Nsplits,4);
corrPredCell = cell(Nsplits,4);
predMatrix   = zeros(height(HSdata),4);

plotTrain  = false;
plotValCV  = false;
plotValTot = false;
plotValTotTot = true;
close all
for i=1:Nsplits
    if plotValCV
        figure
    end
    
    for aa=1:4
        % *** estimate optimal class threshold on training set ***
        activ = CVresults.train.activations{i,aa};
        J     = CVresults.train.J{i,aa};
        % ¤¤ CHOOSE DISEASE AND CLASS THRESHOLD ¤¤
        targetType = "murmur";
        classThr = 2;
        if targetType=="murmur"
            targetVar = sprintf('murGrade%g',aa);
        elseif targetType=="AVMPG"
            targetVar = 'ASPGgrade';
            ASpgThreshold = 10;
            classThr = 1;
            HSdata.ASPGgrade = HSdata.avmeanpg>=ASpgThreshold;
        else
            targetVar = sprintf('%sgrade',targetType);
        end
        Ytarget = HSdata.(targetVar)(J)>=classThr;
        if plotTrain
            subplot(2,2,aa)
        end
        [AUC,X,Y,T] = performanceSummaryNeurNet([],Ytarget,activ,...
                                                [],[],plotTrain);
        AUCmat.train(i,aa) = AUC.whole;
        if plotTrain
            title('training performance curve')
        end
        minSn = 0.92;
        u0 = getOptimalThr(X,Y,T,minSn);
        
        % *** use obtained threshold to estimate SN,SP and AC ***
        activ = CVresults.val.activations{i,aa};
        J     = CVresults.val.J{i,aa};
        Ytarget = HSdata.(targetVar)(J)>=classThr;
        if plotValCV
            subplot(2,2,aa)
        end
        [AUC,X,Y,~,~] = performanceSummaryNeurNet([],Ytarget,activ,...
                                                  [],[],plotValCV);
        
        AUCmat.val(i,aa) = AUC.whole;
        predictions{i,aa} = activ>=u0;
        activations{i,aa} = activ;
        Ytrue{i,aa}       = Ytarget;
        corrPredCell{i,aa} = (activ>=u0)==(Ytarget);
        sn(i,aa) = condProb(activ>=u0,Ytarget);
        sp(i,aa) = condProb(activ<u0,~Ytarget);
        ac(i,aa) = mean((activ>=u0)==Ytarget);
    end
    
    J = cell2mat(CVresults.val.J(i,:)');
    activ = cell2mat(CVresults.val.activations(i,:)');
    Ytarget = cell2mat(Ytrue(i,:)');
    AUCmat.val(i,5) = getAUC(Ytarget,activ);    
end
perfMeasuresOld = [mean(sn);mean(sp);mean(ac)];
[median(sn)',median(sp)',median(ac)']

if targetType=="murmur"
    cols = 1:5;
else
    cols = 1:4;
end
AUCmat.val = AUCmat.val(:,cols);

% *** plot ROC-curve for all predictions for each location ***
if plotValTot
for aa=1:4
    J       = cell2mat(CVresults.val.J(:,aa))
    activ   = cell2mat(CVresults.val.activations(:,aa));
    Ytarget = HSdata.(targetVar)(J);
    AUC = performanceSummaryNeurNet([],Ytarget,activ,...
                                    [],[],true);
end
end

% *** plot ROC-curve for predictions, all positions combined ***
Ytarget = cell2mat(reshape(Ytrue,[Nsplits*4,1]));
activ   = cell2mat(reshape(CVresults.val.activations,[Nsplits*4,1]));
pred = cell2mat(reshape(predictions,[Nsplits*4,1]));

if classThr==2
    figure
else
    hold on
end
% figure
AUC = performanceSummaryNeurNet([],Ytarget,activ,...
                                [],[],plotValTotTot);
predAccuracyTot = balancedPerfEstUsingCVoutput(pred,Ytarget,1);
% title(sprintf('AUC=%g',AUC.whole*100))
% legend(sprintf('murmur-cutoff=1 (AUC=%g)',AUC.whole), 'murmur-cutoff=2 (AUC=0.969)')
legend(sprintf('murmur-cutoff=1 (AUC=%g)',round(AUC.whole,2)))

% *** estimate AUC and prediction accuracy metrics ***
clear T
AUCci = computeCImeanEst(AUCmat.val,"2");
if plotValTot
    figure
end

for aa=1:4
    % get AUC:
    J = cell2mat(CVresults.val.J(:,aa));
    activ = cell2mat(CVresults.val.activations(:,aa));
    if targetType=="murmur"
        Ytarget = HSdata.(sprintf('murGrade%g',aa))(J)>=classThr;
    else
        Ytarget = HSdata.(targetVar)(J)>=classThr;
    end
    if plotValTot
        subplot(2,2,aa)
    end
    AUC = performanceSummaryNeurNet([],Ytarget,activ,...
                                    [],[],plotValTot);
    
    pred  = cell2mat(predictions(:,aa));
    AUCci = computeCImeanEst(AUCmat.val,"2");
    perf = AUCci(aa,:);
    aucTot = [perf(1);perf(2);perf(3)];
    
    predMetrics = balancedPerfEstUsingCVoutput(pred,Ytarget,1);
    % store in structure:
    thrStr = num2name(classThr)
    accurMatrix = [predMetrics{1}, AUCci(aa,:)'];
    accurMatrix(end+1,:) = (accurMatrix(3,:)-accurMatrix(2,:))/2;
    T{aa} = array2table(round(100*accurMatrix,1),...
                                      'V',{'sn','sp','acc','auc'},...
                                      'R',{'Estimate','ci lower','ci upper','half ci width'});
    T{aa} = giveTitle2table(T{aa},...
            sprintf("MGPA prediction of %s>=%g, AA=%g",targetType,classThr,aa)); %#ok<*SAGROW>
    
end

predPerf.murPred.eachAA.(targetType).(thrStr).AUCmat = AUCmat.val;

% formula for predPerf: predPerf.(net).(predictorVar).(targetVarType).(thrStr){aa}
predPerf.murPred.eachAA.(targetType).(thrStr).T = T;
for aa=1:4
    predPerf.murPred.eachAA.(targetType).(thrStr).T{aa}
end

%% generate JOINT ACTIVATIONS splits by borrowing from training set
new = createJointValAndTrainSets(CVresults);
CVresultsJoint = getJointActivationCells(new);
%% estimate PREDICTIVE power of MAX-MURMUR-PREDICTION algorithm
% ¤¤ SET NUMBER OF CV-SPLITS ¤¤
Nsplits = 8;
% get vector of indeces of all positions with validation predictions:
JjointTot = CVresultsJoint.val{1}(:,1);
if Nsplits>1
    for i=2:8
        JjointTot = union(JjointTot,CVresultsJoint.val{i}(:,1));
    end
end

ActMatVal = CVresults.val.activMat;

% define storage variables:
sn    = zeros(Nsplits,1);
sp    = zeros(Nsplits,1);
ac    = zeros(Nsplits,1);
AUCnew = zeros(Nsplits,1);
predictions = cell(Nsplits,1);
activations = cell(Nsplits,1);
Ytrue = cell(Nsplits,1);
Jall = cell2mat(CVresults.valTot.J);

% allVHDpredMatrix = zeros(height(HSdata),4);
% plot settings:
plotTrain = false;
plotVal   = false;
plotValTot = true;
close all
for i=1:Nsplits
    J     = CVresultsJoint.train{i}(:,1);
    activ = CVresultsJoint.train{i}(:,2);
    % set which type of prediction based on murmur grade to use:
    estimatePredictiveModel = true;
    AUCmaximization         = false;
    targetType = "AS";
    classThr      = 1;
    if targetType=="murmur"
        targetVar = 'maxMeanMurGrade';
    else
        targetVar = sprintf('%sgrade',targetType);
    end
    Ytarget = HSdata.(targetVar)(J)>=classThr;
    if plotTrain
        figure
    end
    ActMatVal = getZeroPaddedActivMatrix(CVresults.val.activations(i,:),...
                                   CVresults.val.J(i,:),height(HSdata));
    ActMatTrain = getZeroPaddedActivMatrix(CVresults.train.activations(i,:),...
                                   CVresults.train.J(i,:),height(HSdata));
    ActMat = ActMatTrain + ActMatVal; % can do this since no overlap
    
    minSn = 0.10;
    minSp = 0.10;
    if estimatePredictiveModel && targetType=="AS"
        plotVal = true;
        [activ,u0,Ytarget,AUC] = get_ASactivations(ActMat,HSdata,...
                                 CVresults.trainTot.I{i},CVresults.valTot.I{i},...
                                    classThr,plotVal,minSn,minSp);
    else
        activ = max(ActMat(J),[],2);
        Ytarget  = HSdata.(targetVar)(J)>=classThr;
    end
    
    [~,X,Y,T,~] = performanceSummaryNeurNet([],Ytarget,activ,...
                                            [],[],plotTrain);
    minSn = 0.2;
    minSp = 0.5;
    u0 = getOptimalThr(X,Y,T,minSn,minSp);

    % *** estimate sensitivity and specificity ***
    J = CVresultsJoint.val{i}(:,1);
    
    Ytarget  = HSdata.(targetVar)(J)>=classThr;
    if plotVal
        figure
    end
    [AUC,X,Y] = performanceSummaryNeurNet([],Ytarget,activ,...
                                          [],[],plotVal); %#ok<*ASGLU>
    AUCnew(i)      = AUC.whole;
    predictions{i} = activ>=u0;
    activations{i} = activ;
    Ytrue{i}       = Ytarget;
    
    ac(i) = mean((activ>=u0)==Ytarget);
    sn(i)  = condProb(activ>=u0,Ytarget);
    sp(i)  = condProb(activ<u0,~Ytarget);
end

predictions = cell2mat(predictions);
activations = cell2mat(activations);
YtargetTot = HSdata.(targetVar)(Jall)>=classThr;

if targetType=="AR"
    allVHDpredMatrix(Jall,1) = predictions;
elseif targetType=="MR"
    allVHDpredMatrix(Jall,2) = predictions;
elseif targetType=="AS"
    allVHDpredMatrix(Jall,3) = predictions;
elseif targetType=="MS"
    allVHDpredMatrix(Jall,4) = predictions;
end

% Ytrue = HSdata.(targetVar)(JjointTot)>=classThr;
if plotValTot
    hold on
    performanceSummaryNeurNet([],YtargetTot,activations,...
                              [],[],plotValTot);
end

% *** investigate missed cases ***
IfalseNeg = and(predictions==0,YtargetTot==1);
IfalsePos = and(predictions==1,YtargetTot==0);
ItruePos = and(predictions==1,YtargetTot==1);
JtruePos = Jall(ItruePos);
JfalseNeg = Jall(IfalseNeg);
JfalsePos  = Jall(IfalsePos);
VHDgradeFalseNeg = HSdata(JfalseNeg,:).(targetVar);
murmurGradeFalseNeg = HSdata(JfalseNeg,:).maxMeanMurGrade;
if targetVar=="ASgrade"
    AVMPGfalseNeg = HSdata(JfalseNeg,:).AVMEANPG_T72';
    AVMPGfalsePos = HSdata(JfalsePos,:).AVMEANPG_T72';
end


% *** estimate AUC and prediction accuracy metrics ***
% get AUC:
AUCci   = computeCImeanEst(AUCnew,"tdist");
aucTot = [mean(AUCnew);AUCci(1);AUCci(2)];
% get SN,SP, and AC:
predMetrics = balancedPerfEstUsingCVoutput(predictions,YtargetTot,1);
% store in structure:
clear T
thrStr = num2name(classThr);
AUCstruct.(net).maxMur.(targetType).(thrStr) = AUC;
T = array2table([predMetrics{1},aucTot],'v',{'sn','sp','acc','auc'})
T = giveTitle2table(T,...
                    sprintf("MMPA prediction of %s>=%g",targetType,classThr))
predPerf.(net).maxMur.(targetType).(thrStr) = T;
predPerf.(net).maxMur.(targetType).(thrStr)



% computeCImeanEst(AUCnew,"2")
%% COMPARE MAX vs INDIVIDUAL PREDICTIONS:
targetVar = 'ASGRADE_T72';
% predVar = 'maxMeanMurGrade';
meanPredPerfMaxPred.(targetVar)
meanPredPerfAllPos.(targetVar)
%% performance of MULTIVARIABLE model that uses MAX-MURMUR-PREDICTIONS
%% preliminary:
% here I want to estimate the performance of the model which takes into
% account several different HEALTH variables.
new = createJointValAndTrainSets(CVresults);
CVresultsJoint = getJointActivationCells(new);
CVjointValMat = cell2mat(CVresultsJoint.val);
JvalJoint = CVjointValMat(:,1);

Nsplits = 8;
% ¤¤ SELECT MURMUR VARIABLE ¤¤
% murVar = 'maxMeanMurGrade';
murVar = 'predMaxMurGrade';
names.all.var = {murVar,'AGE_T7','PULSESPIRO_T72','SEX_T7',...
                'DYSPNEA_FAST_UPHILL_T7','CHEST_PAIN_FAST_T7','HIGH_BLOOD_PRESSURE_T7',...
                'DIABETES_T7','SMOKE_DAILY_Q2_T7'};
names.all.categorical = {'SEX_T7','DYSPNEA_FAST_UPHILL_T7','CHEST_PAIN_FAST_T7',...
                'HIGH_BLOOD_PRESSURE_T7','DIABETES_T7','SMOKE_DAILY_Q2_T7'};
names.AR.var = {murVar,'AGE_T7','SEX_T7','DYSPNEA_FAST_UPHILL_T7'};
names.MR.var = {murVar,'AGE_T7','SEX_T7','DYSPNEA_FAST_UPHILL_T7','PULSESPIRO_T72'};
names.AS.var = {murVar,'AGE_T7','SEX_T7'};
names.MS.var = {murVar,'AGE_T7','PULSESPIRO_T72'};

names.AR.model = getLinearModelFormula(names.AR.var,'AR');
names.MR.model = getLinearModelFormula(names.MR.var,'MR');
names.AS.model = getLinearModelFormula(names.AS.var,'AS');
names.MS.model = getLinearModelFormula(names.MS.var,'MS');

HSdata = HSdataTrain;
for i=4:numel(names.all.var)
    HSdata.(names.all.var{i}) = categorical(HSdata.(names.all.var{i}));
end

disease  = 'AS';
classThr = 1;
analyzeErrors = true;

% storage variables:
targetVar  = sprintf('%sgrade',disease);
activations    = cell(Nsplits,1);
murPredictions = cell(Nsplits,1);
predictions    = cell(Nsplits,1);
Jusable        = cell(Nsplits,1);
Ytrue          = cell(Nsplits,1);

close all
plotVal   = false;
plotTrain = false;
plotTot   = true;
if plotVal || plotTrain
    figure
end
for i=1:Nsplits
    % *** training ***
    Jtrain = CVresultsJoint.train{i}(:,1);
    if murVar=='predMaxMurGrade'  
        murActiv = CVresultsJoint.train{i}(:,2);
        dataTrain = [HSdata(Jtrain,:), table(murActiv,'v',{murVar})];
    else
        murActiv  = HSdata(Jtrain,:).(murVar);
        dataTrain = HSdata(Jtrain,:);
    end
    
    % remove rows with nan values:
    Inan = (ones(1,height(dataTrain))==0);
    for j=1:numel(names.(disease).var)
        if ismember(names.(disease).var{j},names.all.categorical)
            Inan = isundefined(dataTrain.(names.(disease).var{j}));
        end
    end
    Inan = (Inan>=1);

    HS    = dataTrain(~Inan,:);
    sigDisease = table(HS.(targetVar)>=classThr,'var',{disease});
    HS    = [HS, sigDisease]; %#ok<*AGROW>
    
    weights = ones(height(HS),1).*(1 + 0*HS.ASGRADE_T72);
    formula = names.(disease).model;
    glm     = fitglm(HS,formula,'distr','binomial',...
                                'Weights',weights);
    
    predTrain  = predict(glm,HS);
    YtrueTrain = HS.(disease);
    [X,Y,T,AUC] = perfcurve(YtrueTrain, predTrain, true);
    if plotTrain
        subplot(2,4,i)
        plot(X,Y)   
    end
    
    % *** validation ***
    Jval = CVresultsJoint.val{i}(:,1);
    if murVar=='predMaxMurGrade' 
        murActiv = CVresultsJoint.val{i}(:,2);
        % append murmur prediction data:
        dataVal  = [HSdata(Jval,:), table(murActiv,'v',{murVar})];
    else
        murActiv = HSdata.(murVar);
        dataVal  = HSdata(Jval,:);
    end
    
    % remove rows with nan values:
    Inan = (ones(1,height(dataVal))==0);
    for j=1:numel(names.(disease).var)
        variable = names.(disease).var{j};
        
        if iscategorical(dataVal.(variable))
            Inan = ismissing(dataVal.(variable));
        else
            Inan = isnan(dataVal.(variable));
        end
    end
    Inan    = (Inan>=1);
    HS = dataVal(~Inan,:);
    
    % append positive class index vector to data frame:
    sigDisease = table(HS.(targetVar)>=classThr,'var',{disease});
    HS    = [HS, sigDisease];
    
    predVal  = predict(glm,HS);
    YtrueVal = HS.(disease);
    if sum(HS.(disease))==0
        AUC = nan;
    else
        [X,Y,T,AUC] = perfcurve(YtrueVal, predVal, true);
        if plotVal
            subplot(2,4,i)
            plot(X,Y)
        end
    end
    
    minSens = 0;
    [perf,predictions{i}] = getPredMetricsAndPred(predTrain,predVal,...
                                               YtrueTrain,YtrueVal,minSens);
    
    % save variables:
    Jusable{i}  = find(~Inan);
    Ytrue{i} = HS.(disease);
    murPredictions{i} = murActiv;
    activations{i} = predVal;
    multiVarValPerf.auc(i) = AUC;
end

% *** performance for all validation set predictions combined ***
valMat      = cell2mat(CVresultsJoint.val);
activations = cell2mat(activations);
murPredictions = cell2mat(murPredictions);
predictions = cell2mat(predictions);
Ytrue       = cell2mat(Ytrue);
J           = valMat(:,1);
if plotTot
    figure
    performanceSummaryNeurNet([],Ytrue,activations,...
                                 numel(activations),[],plotTot);
end

if analyzeErrors
    % locate false negatives (missed positive cases):
    Jerr = find(and(Ytrue==0,predictions==1));
    IndxMissedCases     = JvalJoint(Jerr);
    VHDgradeMissedCases = HSdata.ASGRADE_T72(IndxMissedCases);
    MGmissedCases       = HSdata.maxMeanMurGrade(IndxMissedCases);
    predMGmissedCases   = murPredictions(Jerr);
    missedCases = table(VHDgradeMissedCases,MGmissedCases,...
                            predMGmissedCases,IndxMissedCases)
    % j=1293 (HSdata) is the biggest outlier in terms of missmatch
    % between AS severity and murmur grade {MG=0.5, AS=2}.
end

% compute ci:s for AUC:
[ciAUC,AUC] = balancedPerfEstUsingCVoutput(activations,Ytrue,Nsplits);

% compute ci:s for predicition metrics:
predMetrics = balancedPerfEstUsingCVoutput(predictions,Ytrue,1)

perfMatrix = [predMetrics{1},[ciAUC.estimate;ciAUC.ciLower;ciAUC.ciUpper]];
perfTable = array2table(perfMatrix,'v',{'sn','sp','ac','auc'},'r',...
                            {'estimate','ci lower','ci upper'});
perfTable = giveTitle2table(perfTable,...
                sprintf("multivariate prediction of %s>=%g",disease,classThr))

                        
if murVar=="predMaxMurGrade"
    predictor = "predMaxMurMultiVar";
else
    predictor = "maxMurMultiVar";
end
clear T
thrStr = num2name(classThr);
T = perfTable;
predPerf.(net).(predictor).(disease).(thrStr) = T;
predPerf.(net).(predictor).(disease).(thrStr)
AUCstruct.(net).(predictor).(disease).(thrStr) = AUC;

predPerf.(net).(predictor).(disease).(thrStr)
missedCases;
%% compare performance maxMurmur predictions with P-VALUES
predictor = "maxMur"
disease = "murmur";
thrStr  = "g3";
pval = pValue(AUCstruct.reg.(predictor).(disease).(thrStr) - ...
              AUCstruct.G2.(predictor).(disease).(thrStr));

predPerf.pValues.(predictor).(disease).(thrStr) = pval
predPerf.pValues.(predictor).(disease).(thrStr)
%% compare performance all locations
predictor = "eachAA"
disease   = "AS";
thrStr    = "g2";
aa = 4;
pval = pValue(AUCstruct.reg.(predictor).(disease).(thrStr){aa} - ...
              AUCstruct.G2.(predictor).(disease).(thrStr){aa});

predPerf.pValues.(predictor).(disease).(thrStr){aa} = pval
predPerf.pValues.(predictor).(disease).(thrStr){aa}


