%% script structure
% * train network on the entire training/developement set, setting aside a
% validation set for training stoppage.

% * obtaining test-set results.

%% Train network, or get training or validation set
%#ok<*NOPTS>
%#ok<*DISPLAYPROG>
load nodesAll
nodes0 = nodesAll;

% *** define some storage variables ***
CVresults.val.I           = cell(1,4);
CVresults.val.activations = cell(1,4);
CVresults.valTot.I        = cell(1,1);
CVresults.train.I           = cell(1,4);
CVresults.train.activations = cell(1,4);
CVresults.trainTot.I        = cell(1,1);
CVresults.develop.I           = cell(1,4);
CVresults.develop.activations = cell(1,4);
CVresults.developTot.I        = cell(1,1);
CVresults.test.I           = cell(1,4);
CVresults.test.activations = cell(1,4);
CVresults.testTot.I        = cell(1,1);


Xtrain = cell(1,4);
Ytrain = cell(1,4);
Xval = cell(1,4);
Yval = cell(1,4);
Xtest = cell(1,4);
Ytest = cell(1,4);

% ¤¤ CHOOSE IF TRAIN NETWORK ¤¤
trainNet = true;
% ¤¤ CHOOSE IF GET TRAINING DATA ¤¤
getTrainData = true;
rng(1)
%% get data, collect info, and train net
% *** get training and test set indeces ***
Itrain0 = findInd(1:height(HSdata),Jtrain0);
Ival0   = findInd(1:height(HSdata),Jval0);
Itest0  = findInd(1:height(HSdata),Jtest0);
for aa=1:4
    display(aa) 
    % remove rows corresponding to noissy observations:
    Iclean   = HSdata.(sprintf('MURMUR_%gNOISE_REF_T72',aa))==0;
    
    Itrain = and(Itrain0,Iclean);
    Ival   = and(Ival0,Iclean);
    Itest  = and(Itest0,Iclean);
    Idevelop = or(Itrain,Ival);
    Jtrain = find(Itrain);
    Jval   = find(Ival);
    Jtest  = find(Itest);
    Jdevelop = find(Idevelop);
    
    % segment extraction settings:
    N.ds = 20; % how much to downsample signal
    N.cs = 4;  % How many heartbeats per segment
    N.os = 2;  % by how many heartbeats do the segments overlap
    N.segPerPCG = 6; % how many segments to extract per recording
    Ncomp = [13,200]; % MFCC of segments are reshaped to these dimensions
    
    % define the labels:
    targetStr = sprintf('murGrade%g',aa);
    Y0    = HSdata.(targetStr);
    
    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
    
    % *** set balance settings ***
    posThr = 1;
    if trainNet
        balanceTrain = true;
    else
        balanceTrain = false;
    end
    balanceVal = true;
    
    if trainNet || getTrainData
       % *** get training data ***
       [Xtrain{aa},Ytrain{aa}] = genTrainOrValSet_new(HSdata,Y0,Jtrain,N,aa,...
                                        balanceTrain,nodes0,Ncomp,[],posThr);
    end
   % *** get validation data ***
   [Xval{aa},Yval{aa}] = genTrainOrValSet_new(HSdata,Y0,Jval,N,aa,...
                                    balanceVal,nodes0,Ncomp,[],posThr);
    % *** get validation data ***
    balanceTest = false;
   [Xtest{aa},Ytest{aa}] = genTrainOrValSet_new(HSdata,Y0,Jtest,N,aa,...
                                             balanceTest,nodes0,Ncomp,[]);
    
    % save indeces and id's for validation and data:
    CVresults.val.I{1,aa}     = Ival;
    CVresults.test.I{1,aa}    = Itest;
    CVresults.train.I{1,aa}   = Itrain;
    CVresults.develop.I{1,aa} = Idevelop;
    CVresults.val.J{1,aa}     = find(Ival);
    CVresults.test.J{1,aa}    = find(Itest);
    CVresults.train.J{1,aa}   = find(Itrain);
    CVresults.develop.J{1,aa} = find(Idevelop);
    
    NumHSrec.val(aa)   = sum(Ival);
    NumHSrec.test(aa)  = sum(Itest);
    NumHSrec.train(aa) = sum(Itrain);
end

CVresults.valTot.I   = unionIterated(CVresults.val.I(1,:),"logical");
CVresults.testTot.I  = unionIterated(CVresults.test.I(1,:),"logical");
CVresults.trainTot.I = unionIterated(CVresults.train.I(1,:),"logical");
CVresults.developTot.I = unionIterated(CVresults.develop.I(1,:),"logical");

% combine into 1 training set and validation set:
if trainNet
    Xtrain1234 = [Xtrain{1};Xtrain{2};Xtrain{3};Xtrain{4}];
    Ytrain1234 = [Ytrain{1};Ytrain{2};Ytrain{3};Ytrain{4}];
end
Xval1234 = [Xval{1};Xval{2};Xval{3};Xval{4}];
Yval1234 = [Yval{1};Yval{2};Yval{3};Yval{4}];

%% $$$ $$$ $$$-- train network --$$$ $$$ $$$ 
if trainNet
    % *** define network architecture ***
    numFeatures = height(Xtrain{1}{1}); % same as input size...
    inputSize   = numFeatures; % number of timeseries to take as input
    numClasses = 2;
    if islogical(Y0)
        layers = [...
            sequenceInputLayer(numFeatures)
            lstmLayer(50)
            lstmLayer(50,'OutputMode','last')
            dropoutLayer(.5)
            fullyConnectedLayer(30)
            fullyConnectedLayer(numClasses)
            softmaxLayer
            classificationLayer];
    elseif isnumeric(Y0)
        layers = [...
            sequenceInputLayer(inputSize)
            lstmLayer(50)
            lstmLayer(50,'OutputMode','last')
            dropoutLayer(.5)
            fullyConnectedLayer(30)
            reluLayer
            fullyConnectedLayer(1)
            regressionLayer];
    end
    
    % *** specify training options ***
    % ¤¤ CHOOSE IF CHECK VALIDATION ACCURACY DURING TRAINING ¤¤
    checkValAccuracy = true;
    miniBatchSize = 2^5;
    validationFrequency = floor(numel(Xtrain1234)/miniBatchSize);
    if checkValAccuracy
        validationData = {Xval1234,Yval1234};
    else
        validationData = [];
    end
    % ¤¤ CHOOSE TRAINING OPTIONS ¤¤
    trainTime = 270*60;
    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',50, ...
        'MiniBatchSize',miniBatchSize, ...
        'GradientThreshold',1, ...
        'Plots','training-progress',...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.5, ...
        'LearnRateDropPeriod',5, ...
        'initialLearnRate',0.002, ...
        'ValidationData',validationData, ...
        'ValidationFrequency',validationFrequency, ...
        'LearnRateDropPeriod',5, ...
        'OutputFcn',@(info)myCostumOutputFcn(info));
    
    % *** train and save network ***
    net = trainNetwork(Xtrain1234,Ytrain1234,layers,options);
end

%% *** get training activations ***
plotVal = true;
close all
figure
for aa=1:4
    Ntrain = NumHSrec.train(aa);
    Xtrain{aa} = getUnbalancedSet(Xtrain{aa},Ntrain);
    Ytrain{aa} = getUnbalancedSet(Ytrain{aa},Ntrain);
    Ypred = predict(net, Xtrain{aa});
    
%     YpredTarget = Yval{aa}>=2;
    YpredTarget = HSdata.ASgrade(CVresults.train.I{aa})>=1;
    if plotVal
        subplot(2,2,aa)
    end
    [AUC,~,~,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                                         Ntrain,N.segPerPCG,plotVal);
    if plotVal
        title(sprintf('AUC test set = %g, pos. %g',AUC.whole,aa))
    end
    CVresults.train.activations{1,aa} = p.murWhole;
end

%% *** get validation activations ***
plotVal = true;
close all
figure
for aa=1:4
    Nval = NumHSrec.val(aa);
    Xval{aa} = getUnbalancedSet(Xval{aa},Nval);
    Yval{aa} = getUnbalancedSet(Yval{aa},Nval);
    Ypred = predict(net, Xval{aa});
    
%     YpredTarget = Yval{aa}>=2;
    YpredTarget = HSdata.ASgrade(CVresults.val.I{aa})>=1;
    if plotVal
        subplot(2,2,aa)
    end
    [AUC,~,~,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                                         Nval,N.segPerPCG,plotVal);
    if plotVal
        title(sprintf('AUC test set = %g, pos. %g',AUC.whole,aa))
    end
    CVresults.val.activations{1,aa} = p.murWhole;
end


%% Test set
plotTest = true;
close all
figure
for aa=1:4
    Ntest = NumHSrec.test(aa);
    Ypred = predict(net, Xtest{aa});
    YpredTarget = Ytest{aa}>=2;
%     YpredTarget = HSdata.ASgrade(CVresults.test.I{aa})>=1;
    if plotTest
        subplot(2,2,aa)
    end
    [AUC,~,~,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                                         Ntest,N.segPerPCG,plotTest);
    if plotTest
        title(sprintf('AUC test set = %g, pos. %g',AUC.whole,aa))
    end
    CVresults.test.activations{1,aa} = p.murWhole;
end
%% Combined test-set prediction of Murmur grade
targetArray = getTruthCellArray(HSdata, CVresults.test.J, "murmur", 1);
activ   = cell2mat(CVresults.test.activations(:) );
Ytarget = cell2mat(targetArray(:));
close all
figure
AUC1 = performanceSummaryNeurNet([],Ytarget,activ,...
                                         Ntest,N.segPerPCG,plotTest);
hold on
targetArray = getTruthCellArray(HSdata, CVresults.test.J, "murmur", 2);
activ   = cell2mat(CVresults.test.activations(:) );
Ytarget = cell2mat(targetArray(:));
AUC2 = performanceSummaryNeurNet([],Ytarget,activ,...
                                         Ntest,N.segPerPCG,plotTest);
                                     
legend(sprintf('murmur-cutoff=1 (AUC=%g)',round(AUC1.whole,3)*100), ...
    sprintf('murmur-cutoff=2 (AUC=%g)',round(AUC2.whole,3)*100))

%% get developement-set info and activation matrices:
CVresults.train.activMat = getZeroPaddedActivMatrix(CVresults.train.activations,...
                                             CVresults.train.J,height(HSdata));
CVresults.val.activMat = getZeroPaddedActivMatrix(CVresults.val.activations,...
                                             CVresults.val.J,height(HSdata));
CVresults.test.activMat = getZeroPaddedActivMatrix(CVresults.test.activations,...
                                             CVresults.test.J,height(HSdata));
                                         
CVresults.develop.activMat = CVresults.train.activMat + CVresults.val.activMat;

%% predict AS with AVMPG-calibrated model
ActMat = CVresults.develop.activMat + CVresults.test.activMat;
plotTest = true
[activ,u0,Ytarget,AUC] = get_ASactivations(ActMat,HSdata,...
                            CVresults.developTot.I, CVresults.testTot.I, plotTest) %#ok<*ASGLU>
                        
% looking at false negative AS-predictions:
avmpgFalsePos = HSdata.avmeanpg(CVresults.testTot.J(and(activ.val>=u0,~Ytarget.val)));
sum(avmpgFalsePos>=10)
mean(avmpgFalsePos>=10)

%% predict significant disease with max-murmur and AVMPG-calibrated model
ActMat = CVresults.develop.activMat + CVresults.test.activMat;
disease = "AR"
classThr = 3;
plotTest = true;
[activ,u0,Ytarget,Ypred,AUC] = get_sigVHDactivations(ActMat,HSdata,CVresults.developTot.I,...
                                                    CVresults.testTot.I,...
                                                    disease,classThr,plotTest);

% loop over diseases:
diseaseNames = {'AR','MR','AS','MS'};
predMatrix   = zeros(sum(CVresults.testTot.I),4);
targetMatrix = zeros(sum(CVresults.testTot.I),4);

for i=1:4
    disease = diseaseNames{i};
    if or(disease=='AS',disease=='MS')
        classThr = 1;
    else
        classThr = 4;
    end
    [activ,u0,Ytarget,Ypred,AUC] = get_sigVHDactivations(ActMat,HSdata,...
                                                   CVresults.developTot.I,...
                                                   CVresults.testTot.I,...
                                                   disease,classThr);
    predMatrix(:,i)   = Ypred;
    targetMatrix(:,i) = Ytarget;
    
    predPerf.(disease).Ncaught = sum(and(Ypred,Ytarget));
    predPerf.(disease).Nmissed = sum(and(~Ypred,Ytarget));
    predPerf.(disease).sn = condProb(Ypred,Ytarget);
    predPerf.(disease).sn = condProb(~Ypred,~Ytarget);
end


Ypred = sum(predMatrix,2)>0;
Ytarget = sum(targetMatrix,2)>0;

predPerf.all.sn = condProb(Ypred,Ytarget);
predPerf.all.sp = condProb(~Ypred,~Ytarget);

predPerf.AR.N.truePos = sum(and(predMatrix(:,1),targetMatrix(:,1)));
sum(and(predMatrix(:,2),targetMatrix(:,2)))
sum(and(predMatrix(:,3),targetMatrix(:,3)))
sum(and(predMatrix(:,4),targetMatrix(:,4)))

%% predicting AR, MR, and MS
ActMat = CVresults.develop.activMat + CVresults.test.activMat;
Itrain = CVresults.developTot.I;
Ival = CVresults.testTot.I;
predMatrix = zeros(sum(CVresults.developTot.I), 4);

maxActMat = max(ActMat,[],2);
targetVar = "MR";
classThr = 3;
activ = maxActMat(CVresults.developTot.I);
Ytarget = HSdata.(sprintf('%sgrade',targetVar))(CVresults.developTot.I)>=classThr;
u0 = getOptimalThr([],Ytarget,activ);

activ = maxActMat(CVresults.testTot.I);
Ytarget = HSdata.(sprintf('%sgrade',targetVar))(CVresults.testTot.I)>=classThr;
Ypred = activ>=u0;
AUC = getAUC(Ytarget,activ);



testSet.(targetVar).AUC = AUC;
testSet.(targetVar).sens = condProb(Ypred,Ytarget);
testSet.(targetVar).spec = condProb(~Ypred,~Ytarget)
testSet.(targetVar).u0 = u0;
testSet.AR











%% Predict AS (manually):
load CVresults_wholeTrainingSet_valStop.mat

ActMat = CVresults.develop.activMat + CVresults.test.activMat;
Itrain = CVresults.developTot.I;
Ival = CVresults.testTot.I;
targetType = 'sigVHD31';
minSn = .5;
minSp = 0;

plotVal = true;
classThr = 1;
[activ,u0,Ytarget,Ypred,AUC] = get_sigVHDactivations(ActMat,HSdata,...
                                         Itrain,Ival,targetType,classThr,...
                                         plotVal,minSn,minSp);

% [activ,u0,Ytarget,Ypred,AUC] = get_ASactivations(ActMat,HSdata,...
%                                                 Itrain,Ival,classThr,plotVal);

Z = succAndFailureAnalysis(Ypred.val,Ytarget.val,Ival,HSdata,...
                {'avmeanpg','ARgrade','MRgrade','ASgrade','MSgrade'})

%%
Ncaught = sum(and(Ypred.val,Ytarget.val))
sum(and(Ypred.val,Ytarget.val))

getAUC(HSdata.avmeanpg(Ival)>10,activ.val)

