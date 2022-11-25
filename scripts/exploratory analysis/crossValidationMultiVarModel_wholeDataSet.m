% data = HSdata0;
% data = HSdata(union(Jtrain0,Jval0),:);
data = HSdata;
data.SEX_T7                 = categorical(data.SEX_T7);
data.HIGH_BLOOD_PRESSURE_T7 = categorical(data.HIGH_BLOOD_PRESSURE_T7);
data.DYSPNEA_FAST_UPHILL_T7 = categorical(data.DYSPNEA_FAST_UPHILL_T7);
data.CHEST_PAIN_FAST_T7     = categorical(data.CHEST_PAIN_FAST_T7>=1);
data.DIABETES_T7            = categorical(data.DIABETES_T7==1);

murAll = {'MGpos1','MGpos2','MGpos3','MGpos4'};
murMax = {'maxMeanMurGrade'};
names.all.categorical = {'SEX_T7','DYSPNEA_FAST_UPHILL_T7','CHEST_PAIN_FAST_T7',...
                'HIGH_BLOOD_PRESSURE_T7','DIABETES_T7'};
names.AR.var = [murMax,{'AGE_T7','SEX_T7','DYSPNEA_FAST_UPHILL_T7'}];
names.MR.var = [murMax,{'AGE_T7','SEX_T7','DYSPNEA_FAST_UPHILL_T7','PULSESPIRO_T72'}];
names.AS.var = [murMax,{'AGE_T7','SEX_T7'}];
names.MS.var = [murMax,{'AGE_T7','PULSESPIRO_T72'}];
names.AVMPG.var = [murMax,{'AGE_T7','SEX_T7','DIABETES_T7'}];

% ¤¤ SELECT DISEASE TO PREDICT ¤¤
% disease = 'AS';
disease = 'MR';
classThr = 3;
varNames = names.(disease).var;
% ¤¤ SELECT CLASS THRESHOLD ¤¤
if disease=="AVMPG" %#ok<*BDSCI>
    targetVarName = 'AVMAXPG_T72';
else
    targetVarName = sprintf('%sgrade',disease);
end
modelFormula = getLinearModelFormula(varNames,disease);
% remove rows with nan values:
clear nanLoss
Inan = (ones(height(data),1)==0);
for j=1:numel(varNames)
    varName   = names.(disease).var{j};
    IposCases = data.ASGRADE_T72>=1;
    variable  = data.(varName);
    if iscategorical(variable)
        nanLoss.(varName) = sum(and(ismissing(variable),IposCases));
        Inan = Inan + ismissing(variable);
    else
        nanLoss.(varName) = sum(and(ismissing(variable),IposCases));
        Inan = Inan + isnan(variable);
    end
end
Inan = Inan + isnan(data.(targetVarName));
Inan = (Inan>=1);
data.(targetVarName)(Inan)
modelData = data(~Inan,:);

% number of CV-partitions:
K = 8;
genNewCV = true;
cvPartitions = cvpartition(modelData.(targetVarName),'Kfold',K);
categoricalTarget = numel(unique(modelData.(targetVarName)))<6;
if ~categoricalTarget
    Y = modelData.(targetVarName);
else
    Y = modelData.(targetVarName)>=classThr;
end
modelData = [modelData,table(Y,'var',{disease})];
MGtable = table(modelData.Murmur_1_grade_ref_ny_T72,'v',{'MGpos1'});
for aa=2:4
    MGtable = [MGtable, ...
        table(modelData.(sprintf('Murmur_%g_grade_ref_ny_T72',aa)),...
                        'var',{sprintf('MGpos%g',aa)} )]; %#ok<*AGROW>
end
modelData = [modelData,MGtable];

predictions = cell(8,1);
activations = cell(8,1);
groundTruth = cell(8,1);
plotIt     = false;
plotAllVal = true;
close all
for i=1:8
    if genNewCV
        Jtrain = cvPartitions.training(i);
        Jval   = cvPartitions.test(i);
    else
        Jtrain = CVresultsJoint.train{i}(:,1);
        Jval   = CVresultsJoint.val{i}(:,1);
        % *** map indeces so they are wrt. to new data-set ***
        J0     = 1:numel(union(Jtrain0,Jval0));
        Jclean = ~Inan;
        Jtrain = getIndecesWRTnewBasis(J0,Jtrain,Jclean);
        Jval   = getIndecesWRTnewBasis(J0,Jval,Jclean);
    end
    modelDataTrain = modelData(Jtrain,:);
    modelDataVal   = modelData(Jval,:);

    % *** train ***
    if ~categoricalTarget
        weights = ones(numel(Y),1) + 0*Y.^4;
        glm = fitlm(modelDataTrain,modelFormula);
    else
        glm = fitglm(modelDataTrain,modelFormula,'distr','binomial');
    end
    % *** estimate optimal threshold ***
    trainPred  = glm.predict;
    if ~categoricalTarget
        trainTruth = modelDataTrain.ASGRADE_T72>=1;
    else
        trainTruth = modelDataTrain.(sprintf('%sgrade',disease))>=classThr;
    end
    [Xtrain,Ytrain,T,AUCtrain(i)] = perfcurve(trainTruth, trainPred, true);
    u0 = getOptimalThr(Xtrain,Ytrain,T,0);
    % u0 = 0.005;
    % get prediction performance
    valPred  = predict(glm,modelDataVal);
    if ~categoricalTarget
        valTruth = modelDataVal.ASGRADE_T72>0;
    else
        valTruth = modelDataVal.(disease);
    end
%     [AUCval(i),Xval,Yval,~,~] = performanceSummaryNeurNet([],valTruth,valPred,...
%                     numel(activ),[],false);
    if sum(valTruth)==0
        AUCval(i) = nan;
    else
        [Xval,Yval,T,AUCval(i)] = perfcurve(valTruth, valPred, true); %#ok<*SAGROW>
    end
    predictions{i} = valPred>=u0;
    groundTruth{i} = valTruth;
    activations{i} = valPred;

    if plotIt && sum(valPred)>0
        figure
        subplot(2,1,1)
            plot(Xtrain,Ytrain)
        subplot(2,1,2)
            plot(Xval,Yval)
    end
    
end

% *** train on whole training set ***
if ~categoricalTarget
    weights = ones(numel(Y),1).*(1 + 0*Y.^4);
    glm = fitlm(modelData,modelFormula,'Weights',weights);
else
    glm = fitglm(modelData,modelFormula,'distr','binomial');
end

activations = cell2mat(activations);
predictions = cell2mat(predictions);
groundTruth = cell2mat(groundTruth);

SN = condProb(predictions,groundTruth);
SP = condProb(~predictions,~groundTruth);
computeCImeanEst(AUCval,"2")
[ciAUC,AUC] = balancedPerfEstUsingCVoutput(activations,groundTruth,8);

if plotAllVal
    figure
    [Xtot,Ytot,Ttot,AUCtot] = perfcurve(groundTruth, activations, true);
    plot(Xtot,Ytot)
    title(sprintf('AUC=%g, SN=%g, SP=%g',mean(AUC)*100,SN,SP))
end

glm.Coefficients
if ~categoricalTarget
    rms(activations - modelData.AVMAXPG_T72(1:numel(activations)))
end
% sum(Inan)
sum(data.ASGRADE_T72(Inan)>0)
% *** AS ***
% AUC validation set: 95.49  97.07  98.651       (maxMur)
% AUC validation set: 94.47  96.635 98.796       (allMur)
% SN = 86.05; SP: 93.10 (using training-set estimated threshold, allMur)

% *** MS ***
% AUC validation set: (93.126  96.523  99.919) (maxMur)
% SN = 76.92; SP: 93.71 (using training-set estimated threshold)
% SN = 92.31; SP: 90.03 (manually tuning the activation threshold)

% CONCLUSION: The specificity and sensitivity is a bit underwhelming for
% MS, but it is better than it was for the models developed the training
% set. Manually tuning the activation threshold yields a much higher
% sensitivity and specificity. This is likely because there are so few
% cases of MS (14 in total) in each training set, which results in high
% variance in the estimation of the optimal activation threshold and the
% their corresponding performance. There are 46 cases of mild or greater AS

% There are 46 cases of mild or greater AS.

% Including all murmur grades in a multivariable model did NOT produce
% significantly better results than using max-murmur only.
%%

% cases whole modelData set:
MScases.all = sum(HSdata0.MSPRESENCE_T72);
% cases whole data all rows with complete info:
MScases.allCompleteInfo = sum(HSdata0.MSPRESENCE_T72(~Inan));
MScases.nonCorruptAudio = sum(HSdata.MSPRESENCE_T72);
MScases.trainingSet = sum(HSdata.MSPRESENCE_T72(union(Jtrain0,Jval0)))
% 1 (0.77%) of removed rows have MS>=1




