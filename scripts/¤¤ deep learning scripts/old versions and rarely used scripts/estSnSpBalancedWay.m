%#ok<*NOPTS>
%#ok<*NASGU>
%#ok<*STCMP,*BDSCI,*BDSCA>
%% load data
net = "reg";
if net=="reg"
    load CVresults_netMurRegAllPos.mat
elseif net=="G2"
    load CVresults_netMurG2AllPos.mat
elseif net=="net1"
    load CVresults_wholeTrainingSet_net1.mat
elseif net=="net2"
    load CVresults_wholeTrainingSet_net2.mat
end
%% Prediction performance for each individual location
clear AUCmat
% *** estimate sensitivity and specificity for each location ***
% ¤¤ REMEMBER TO SET net = "nameOfNetwork"
if net=="net1" || net=="net2"
    modelData = HSdata;
else
    modelData = HSdata(union(Jtrain0,Jval0),:);
end

% ¤¤ SET NUMBER OF CV-SPLITS ¤¤
Nsplits = 8;
% storage variables:
sn = zeros(Nsplits,4);
sp = zeros(Nsplits,4);
ac = zeros(Nsplits,4); 
AUCmat.val = zeros(Nsplits,4);
AUCmat.train = zeros(Nsplits,4);
predictions  = cell(Nsplits,4);
activations  = cell(Nsplits,4); 
Ytrue        = cell(Nsplits,4);
corrPredCell = cell(Nsplits,4);
predMatrix   = zeros(height(modelData),4);

plotTrain = false;
plotVal   = false;
close all
for i=1:Nsplits
    if plotVal;figure;end
    for aa=1:4
        % *** estimate optimal class threshold on training set ***
        activ = CVresults.train.activations{i,aa};
        J     = CVresults.train.J{i,aa};
        % ¤¤ CHOOSE DISEASE AND CLASS THRESHOLD ¤¤
        targetVarType = "AS";
        classThr = 1;
        if targetVarType=="murmur"
            targetVar = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
        elseif targetVarType=="sigVHD"
            targetVar = targetVarType;
        else
            targetVar = sprintf('%sGRADE_T72',targetVarType);
        end
        Ytarget = modelData.(targetVar)(J)>=classThr;
        if plotTrain
            figure
        end
        [AUC,X,Y,T] = performanceSummaryNeurNet([],Ytarget,activ,...
            numel(activ),[],plotTrain);
        AUCmat.train(i,aa) = AUC.whole;
        if plotTrain
            title('training performance curve')
        end
        minSensitivity = 0.1;
        u0 = getOptimalThr(X,Y,T,minSensitivity);
        
        % *** use obtained threshold to estimate SN,SP and AC ***
        activ = CVresults.val.activations{i,aa};
        J     = CVresults.val.J{i,aa};
        Ytarget = modelData.(targetVar)(J)>=classThr;
        if plotVal
            subplot(2,2,aa)
        end
        [AUC,X,Y,~,~] = performanceSummaryNeurNet([],Ytarget,activ,...
                                                  [],[],plotVal);
        
        AUCmat.val(i,aa) = AUC.whole;
        predictions{i,aa} = activ>=u0;
        activations{i,aa} = activ;
        Ytrue{i,aa}       = Ytarget;
        corrPredCell{i,aa} = (activ>=u0)==(Ytarget);
        sn(i,aa) = condProb(activ>=u0,Ytarget);
        sp(i,aa) = condProb(activ<u0,~Ytarget);
        ac(i,aa) = mean((activ>=u0)==Ytarget);
    end
end
perfMeasuresOld = [mean(sn);mean(sp);mean(ac)];
[median(sn)',median(sp)',median(ac)']

% *** plot ROC-curve for all predictions for each location ***
if plotVal
figure
for aa=1:4
    activ   = cell2mat(CVresults.val.activations(:,aa));
    Ytarget = cell2mat(Ytrue(:,aa));
    subplot(2,2,aa)
    AUC = performanceSummaryNeurNet([],Ytarget,activ,...
            numel(activ),[],true);
end
end

% *** plot ROC-curve for predictions, all positions combined ***
Ytarget = cell2mat(reshape(Ytrue,[Nsplits*4,1]));
activ = cell2mat(reshape(CVresults.val.activations,[Nsplits*4,1]));
hold on
AUC = performanceSummaryNeurNet([],Ytarget,activ,...
        numel(activ),[],true);
    
legend('murmur-cutoff=1 (AUC=0.929)', 'murmur-cutoff=2 (AUC=0.966)')


clear T
for aa=1:4
    % *** estimate AUC and prediction accuracy metrics ***
    % get AUC:
    activ = cell2mat(activations(:,aa));
    truth = cell2mat(Ytrue(:,aa));
    pred  = cell2mat(predictions(:,aa));
    [perf,AUC] = balancedPerfEstUsingCVoutput(activ,truth,Nsplits);
    aucTot = [perf.estimate;perf.ciLower;perf.ciUpper];
    predMetrics = balancedPerfEstUsingCVoutput(pred,truth,1);
    % store in structure:
    thrStr = num2name(classThr)
    AUCstruct.(net).eachAA.(targetVarType).(thrStr){aa} = AUC;
    T{aa} = array2table([predMetrics{1},aucTot],...
                                      'v',{'sn','sp','acc','auc'});
    T{aa} = giveTitle2table(T{aa},...
                  sprintf("MGPA prediction of %s>=%g, AA=%g",targetVarType,classThr,aa)); %#ok<*SAGROW>
    
end

% formula for predPerf: predPerf.(net).(predictorVar).(targetVarType).(thrStr){aa}
predPerf.(net).eachAA.(targetVarType).(thrStr) = T;
for aa=1:4
    predPerf.(net).eachAA.(targetVarType).(thrStr){aa}
end
%% repeat above with MODIFIED CV sets to estimate BIAS introduced
new = createJointValAndTrainSets(CVresults);    

% *** preliminary ***
sn  = zeros(Nsplits,4);
sp  = zeros(Nsplits,4);
ac  = zeros(Nsplits,4);
auc = zeros(Nsplits,4);
Ytrue       = cell(Nsplits,4);
activations = cell(Nsplits,4);
predictions = cell(Nsplits,4);

plotIt = false;
% close all
for i=1:Nsplits
    if plotIt;figure;end
    for aa=1:4
        % *** estimate activation threshold ***
        activ = new.trainMat{i,aa}(:,2);
        J     = new.trainMat{i,aa}(:,1);
        targetVar = 'ASGRADE_T72';
%         predVar = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
        classThr = 1;
        Ytarget = modelData.(targetVar)(J)>=classThr;
        if plotIt
            subplot(4,2,aa)
        end
        [~,X,Y,T,~] = performanceSummaryNeurNet([],Ytarget,activ,...
                                           numel(activ),[],plotIt);
        minSensitivity = 0.9;
        u0 = getOptimalThr(X,Y,T,minSensitivity);

        % *** estimate sensitivity and specificity ***
        activ = new.valMat{i,aa}(:,2);
        J     = new.valMat{i,aa}(:,1);
        Ytarget = modelData.(targetVar)(J)>=classThr;
        if plotIt
            subplot(4,2,4+aa)
        end
        AUC = performanceSummaryNeurNet([],Ytarget,activ,...
                                        numel(activ),[],plotIt);

        Ytrue{i,aa}       = Ytarget;
        activations{i,aa} = activ;
        predictions{i,aa} = activ>=u0;
        ac(i,aa) = mean((activ>=u0)==Ytarget);
        sn(i,aa)  = condProb(activ>=u0,Ytarget);
        sp(i,aa)  = condProb(activ<u0,~Ytarget);
        if isstruct(AUC)
            auc(i,aa) = AUC.whole;
        else
            auc(i,aa) = nan;
        end
    end
end

aucTot = zeros(4,1);
figure
snTot = zeros(4,1);
spTot = zeros(4,1);
accTot = zeros(4,1);
for aa=1:4
    valset = cell2mat(new.valMat(:,aa));
    activ = valset(:,2);
    J     = valset(:,1);
    Ytarget = modelData.(targetVar)(J)>=classThr;
    subplot(2,2,aa)
    AUC = performanceSummaryNeurNet([],Ytarget,activ,...
                                    numel(activ),[],true);
    snTot(aa) = condProb(cell2mat(predictions(:,aa)),Ytarget);
    spTot(aa) = condProb(~cell2mat(predictions(:,aa)),~Ytarget);
    accTot(aa) = prob(cell2mat(predictions(:,aa))==Ytarget);
    aucTot(aa) = AUC.whole;
end

% snTot  = mean(sn,'omitnan')';
% spTot  = mean(sp,'omitnan')';
% accTot = mean(acc,'omitnan')';
meanPredPerfAllPosMod.(targetVar) = table(snTot,spTot,accTot,aucTot,...
                                'RowNames',{'pos1','pos2','pos3','pos4'});
meanPredPerfAllPosMod.(targetVar)
% difference between old and new CV-performance:      sn    sp    acc
% rms(perfNew-perfOld) = rms([0    0.05 0.17 -0.12]) = 0.1040%
% the results are very near identical, and I conclude that bias is
% insignificant.
%% investigate PERFORMANCE DIFFERENCE to see if BIAS is significant
T1 = meanPredPerfAllPos.(targetVar);
T2 = meanPredPerfAllPosMod.(targetVar);
T = table(T1,T2) 
% idSickAndNoisy = find(and(modelData.ASPRESENCE_T72,modelData.MURMUR_2NOISE_REF_T72))

mean(meanPredPerfAllPosMod.(targetVar).aucTot - meanPredPerfAllPos.(targetVar).aucTot)
mean(meanPredPerfAllPosMod.(targetVar).accTot - meanPredPerfAllPos.(targetVar).accTot)
% conclusions:
% AS>0: mean difference for AUC is -0.0015, and -0.0064 for sensitivity,
% and 0.0048 for accuracy.
%% test
aa = 1;
i  = 6;
[new.valMat{i,1}(end,1),new.valMat{i,2}(end,1),...
 new.valMat{i,3}(end,1),new.valMat{i,4}(end,1)]
[new.valMat{i+1,1}(1,1),new.valMat{i+1,2}(1,1),...
 new.valMat{i+1,2}(1,1),new.valMat{i+1,3}(1,1)]
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
% define storage variables:
sn     = zeros(Nsplits,1);
sp     = zeros(Nsplits,1);
ac    = zeros(Nsplits,1);
AUCnew = zeros(Nsplits,1);
u0vec = cell(Nsplits,1);
predictions = cell(Nsplits,1);
activations = cell(Nsplits,1);
% plot settings:
plotTrain = false;
plotVal   = true;
close all
for i=1:Nsplits
    J     = CVresultsJoint.train{i}(:,1);
    activ = CVresultsJoint.train{i}(:,2);
    targetVarType = "AS";
    classThr      = 1;
    if targetVarType=="murmur"
        targetVar = 'maxMeanMurGrade';
    elseif targetVarType=="sigVHD"
        targetVar = targetVarType;
    else
        targetVar = sprintf('%sGRADE_T72',targetVarType);
    end
    Ytarget = modelData.(targetVar)(J)>=classThr;
    if plotTrain
        figure
    end
    [~,X,Y,T,~] = performanceSummaryNeurNet([],Ytarget,activ,...
                numel(activ),[],plotTrain);
    minSensitivity = 0.1;
    u0 = getOptimalThr(X,Y,T,minSensitivity,0.85);
    
    % *** estimate sensitivity and specificity ***
    J     = CVresultsJoint.val{i}(:,1);
    activ = CVresultsJoint.val{i}(:,2);
%     predVar  = 'ARGRADE_T72';
%     classThr = 3;
    Ytarget  = modelData.(targetVar)(J)>=classThr;
    if plotVal
        figure
    end
    [AUC,X,Y] = performanceSummaryNeurNet([],Ytarget,activ,...
                numel(activ),[],plotVal); %#ok<*ASGLU>
    predictions{i} = activ>=u0;
    activations{i} = activ;
    u0vec{i}       = ones(numel(activ),1)*u0;
    ac(i) = mean((activ>=u0)==Ytarget);
    sn(i)  = condProb(activ>=u0,Ytarget);
    sp(i)  = condProb(activ<u0,~Ytarget);
end

u0vec = cell2mat(u0vec);
predictions = cell2mat(predictions);
activations = cell2mat(activations);
predictions = activations>=u0vec;
Ytrue = modelData.(targetVar)(JjointTot)>=classThr;
if plotVal
    figure
    AUC = performanceSummaryNeurNet([],Ytrue,activations,...
                                    [],[],plotVal);
end

% *** investigate missed cases ***
ImissedCases = and(predictions==0,Ytrue==1);
IndxMissedCases = JjointTot(ImissedCases);
modelData(IndxMissedCases,:).ASGRADE_T72
modelData(IndxMissedCases,:).maxMeanMurGrade
% *** estimate AUC and prediction accuracy metrics ***
% get AUC:
[perf,AUC]   = balancedPerfEstUsingCVoutput(activations,Ytrue,Nsplits);
aucTot = [perf.estimate;perf.ciLower;perf.ciUpper];
% get SN,SP, and AC:
predMetrics = balancedPerfEstUsingCVoutput(predictions,Ytrue,1);
% store in structure:
clear T
thrStr = num2name(classThr);
AUCstruct.(net).maxMur.(targetVarType).(thrStr) = AUC;
T = array2table([predMetrics{1},aucTot],'v',{'sn','sp','acc','auc'})
T = giveTitle2table(T,...
                    sprintf("MMPA prediction of %s>=%g",targetVarType,classThr))
predPerf.(net).maxMur.(targetVarType).(thrStr) = T;
predPerf.(net).maxMur.(targetVarType).(thrStr)
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

modelData = HSdataTrain;
for i=4:numel(names.all.var)
    modelData.(names.all.var{i}) = categorical(modelData.(names.all.var{i}));
end

disease  = 'AR';
classThr = 3;
analyzeErrors = true;

% storage variables:
targetVar  = sprintf('%sGRADE_T72',disease);
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
%     Jtrain = CVresultsBalanced.train{i};
    if murVar=='predMaxMurGrade'  
        murActiv  = CVresultsJoint.train{i}(:,2);
        dataTrain = [modelData(Jtrain,:), table(murActiv,'v',{murVar})];
    else
        murActiv  = modelData(Jtrain,:).(murVar);
        dataTrain = modelData(Jtrain,:);
    end
    
    % remove rows with nan values:
    Inan = (ones(1,height(dataTrain))==0);
    for j=1:numel(names.(disease).var)
        if ismember(names.(disease).var{j},names.all.categorical)
            Inan = isundefined(dataTrain.(names.(disease).var{j}));
        end
    end
    Inan = (Inan>=1);

    dataMdl    = dataTrain(~Inan,:);
    sigDisease = table(dataMdl.(targetVar)>=classThr,'var',{disease});
    dataMdl    = [dataMdl, sigDisease]; %#ok<*AGROW>
    
    weights = ones(height(dataMdl),1);
    formula = names.(disease).model;
    glm     = fitglm(dataMdl,formula,'distr','binomial',...
                                'Weights',weights);
    
    predTrain  = predict(glm,dataMdl);
    YtrueTrain = dataMdl.(disease);
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
        dataVal  = [modelData(Jval,:), table(murActiv,'v',{murVar})];
    else
        murActiv = modelData.(murVar);
        dataVal  = modelData(Jval,:);
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
    dataMdl = dataVal(~Inan,:);
    
    % append positive class index vector to data frame:
    sigDisease = table(dataMdl.(targetVar)>=classThr,'var',{disease});
    dataMdl    = [dataMdl, sigDisease];
    
    predVal  = predict(glm,dataMdl);
    YtrueVal = dataMdl.(disease);
    if sum(dataMdl.(disease))==0
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
    Ytrue{i} = dataMdl.(disease);
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
    VHDgradeMissedCases = modelData.ASGRADE_T72(IndxMissedCases);
    MGmissedCases       = modelData.maxMeanMurGrade(IndxMissedCases);
    predMGmissedCases   = murPredictions(Jerr);
    missedCases = table(VHDgradeMissedCases,MGmissedCases,...
                            predMGmissedCases,IndxMissedCases)
    % j=1293 (modelData) is the biggest outlier in terms of missmatch
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


