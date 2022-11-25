% This script trains networks (any target) on the cross validation
% partitions, using data from all four positions. It allows setting either
% a binary or continuous training target.

% Load data on CV-partitioning and segmentation:
load nodesAll.mat
load CVpartitions.mat
nodes0 = nodesAll;

Nsplits = 2;

% ¤¤ CHOOSE WHETHER TO TRAIN NET OR LOAD PRETRAINED NETWORK FOR TESTING ¤¤
trainNet = true;
if trainNet
    % ¤¤ IF TRAIN: CHOOSE IF SAVE NETWORKS ¤¤
    loadNetsName = "";
elseif not(trainNet)
    % ¤¤ IF NOT TRAIN: SELECT NETS TO LOAD ¤¤
    loadNetsName = "networksCVnoNoiseDrOutMurRegAllPosValStopOvertrain.mat"
    load(loadNetsName)
end

% ¤¤ CHOOSE IF GET TRAINING DATA ¤¤
getTrainData = true;
% ¤¤ CHOOSE TYPE OF TARGET THAT NETWORK PREDICTS (murmur,ARgrade,...) ¤¤
targetType = "noise";
% ¤¤ CHOOSE IF RESAMPLE TO BALANCE DATA ¤¤
balanceTrain = true;
balanceVal   = true;
% ¤¤ CHOOSE IF TRAINING OUTPUT VARIABLE IS GRADED OR BINARY ¤¤
gradedOutput = false;
binaryOutput = ~gradedOutput;

if binaryOutput
    % ¤¤ IF BINARY OUTPUT: CHOOSE CLASS SPERATION THRESHOLD ¤¤
    classThr = 1;
    posThr = classThr;
    
elseif gradedOutput
    % ¤¤ IF GRADED OUTPUT: SET THRESHOLD FOR RESAMPLING CLASS ¤¤
    posThr = 1; 
    classThr = 2; % class threshold for plotting
end
% *** note: resampling is done for all Y >= posThr ***
% ¤¤ PLOT SETTINGS ¤¤
plotTrain = false;
plotVal = true;
close all


Info = table(trainNet,getTrainData,balanceTrain,balanceVal,binaryOutput,...
             loadNetsName)


% define some result storage variables:
CVresults.val.I     = cell(Nsplits,4);
CVresults.val.activ = cell(Nsplits,4);
CVresults.train.I     = cell(Nsplits,4);
CVresults.train.activ = cell(Nsplits,4);
CVresults.valTot.I   = cell(Nsplits,1);
CVresults.trainTot.I = cell(Nsplits,1);

% keep track of number of recordings for each validation set and position:
NumHSrec.val   = zeros(4,1);
NumHSrec.train = zeros(4,1);
AUCmat = zeros(Nsplits,5); % the fifth column is for the AUC corresponding to all positions

for i=1:Nsplits
disp(i) %#ok<*NOPTS>
Xtrain = cell(4,1);
Ytrain = cell(4,1);
XvalBal = cell(4,1);
YvalBal = cell(4,1);
Xval   = cell(4,1);
Yval   = cell(4,1);

% get index for rows corresponding to the training and validation sets:
ItrainRows = CVpartitions.train.I{i};
IvalRows   = CVpartitions.val.I{i};

% $$$ $$$ $$$-- get validation and/or training data --$$$ $$$ $$$
for aa=1:4
    disp(aa)
    % get index of clean recordings for position aa:
    if targetType=="noise"
        Iclean = ones(height(HSdata),1);
    else
        Iclean = HSdata.(sprintf('noise%g',aa))==0;
    end
    % get index for training and validation set observations of position aa:
    Itrain = and(ItrainRows,Iclean);
    Ival   = and(IvalRows,Iclean);
    Jtrain = find(Itrain);
    Jval   = find(Ival);
    
    if targetType=="murmur"
        targetStr = sprintf('murGrade%g',aa);
    elseif targetType=="noise"
        targetStr = sprintf('noise%g',aa);
    else
        targetStr = targetType;
    end
        
    if gradedOutput
        Y0 = HSdata.(targetStr);
    elseif binaryOutput
        Y0 = HSdata.(targetStr)>=classThr;
    end
    
    % segment extraction settings:
    N.ds = 20; % how much to downsample signal
    N.cs = 4;  % How many heartbeats per segment
    N.os = 2;  % by how many heartbeats do the segments overlap
    N.segPerPCG = 6;  % how many segments to extract per recording
    Ncomp = [13,200]; % MFCC of segments are reshaped to these dimensions
    
    if trainNet || getTrainData
       % *** get training data ***
       [Xtrain{aa},Ytrain{aa}] = genTrainOrValSet_new(HSdata,Y0,Jtrain,N,aa,...
                                        balanceTrain,nodes0,Ncomp,[],posThr);
    end
    
   % *** get validation data ***
   [XvalBal{aa},YvalBal{aa},Xval{aa},Yval{aa}] = genTrainOrValSet_new(HSdata,Y0,Jval,N,aa,...
                                    balanceVal,nodes0,Ncomp,[],posThr);
   % save indeces and id's for validation and data:
    CVresults.val.I{i,aa}   = Ival;
    CVresults.train.I{i,aa} = Itrain;
    
    NumHSrec.val(aa)   = sum(Ival);
    NumHSrec.train(aa) = sum(Itrain);
end
% save indeces for rows in validation/training set for which there is at
% least one prediction (some rows will have no predictions due to all 4
% recordings being noisy):
CVresults.valTot.I{i,1}   = unionIterated(CVresults.val.I(i,:),"logical");
CVresults.trainTot.I{i,1} = unionIterated(CVresults.train.I(i,:),"logical");

% combine data from all pos. into one training and validation set to feed
% to training algorithm:
Xtrain1234  = [Xtrain{1} ;Xtrain{2} ;Xtrain{3} ;Xtrain{4}];
Ytrain1234  = [Ytrain{1} ;Ytrain{2} ;Ytrain{3} ;Ytrain{4}];
XvalBal1234 = [XvalBal{1};XvalBal{2};XvalBal{3};XvalBal{4}];
YvalBal1234 = [YvalBal{1};YvalBal{2};YvalBal{3};YvalBal{4}];

% $$$ $$$ $$$-- train network --$$$ $$$ $$$ 
if trainNet
    % *** define network architecture ***
    numFeatures = height(Xtrain{1}{1}); % same as input size...
    inputSize   = numFeatures; % number of timeseries to take as input
    numClasses = 2;
    if islogical(Y0)
        layers = [ ...
            sequenceInputLayer(numFeatures)
            lstmLayer(50)
            lstmLayer(50,'OutputMode','last')
            dropoutLayer(.5)
            fullyConnectedLayer(30)
            fullyConnectedLayer(numClasses)
            softmaxLayer
            classificationLayer];
        
    elseif isnumeric(Y0)
        layers = [ ...
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
    checkValAccuracy = false;
    miniBatchSize = 2^5;
    validationFrequency = floor(numel(Xtrain1234)/miniBatchSize);
    if checkValAccuracy
        validationData = {XvalBal1234,YvalBal1234};
    else
        validationData = [];
    end
    % ¤¤ CHOOSE TRAINING OPTIONS ¤¤
    trainTime = 270*60;
    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',25, ...
        'MiniBatchSize',miniBatchSize, ...
        'GradientThreshold',1, ...
        'Plots','training-progress', ...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.5, ...
        'LearnRateDropPeriod',5, ...
        'initialLearnRate',0.002, ...
        'ValidationData',validationData, ...
        'ValidationFrequency',validationFrequency, ...
        'LearnRateDropPeriod',5, ...
        'OutputFcn',@(info)myCostumOutputFcn(info));
    
    % *** train and save network ***
    % ¤¤ CHOOSE IF INITIATE WITH WEIGHTS FROM PRETRAINED NET ¤¤
    initTrainingWithNets = false;
    if initTrainingWithNets
        initNetsName = "networksCVnoNoiseDrOutMurRegAllPosValStop.mat";
        initNets = load(initNetsName);
        initNet = initNets.networks{i};
    end
    InfoTrain = table(checkValAccuracy,initTrainingWithNets)
    
    if initTrainingWithNets
        netTemp = trainNetwork(Xtrain1234,Ytrain1234,initNet.Layers,options);
    else
        netTemp = trainNetwork(Xtrain1234,Ytrain1234,layers,options);
    end
    
    % save nets in folder:
    location = 'C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\matlab files\¤¤ deep learning scripts\¤¤ trained networks\tempFiles';
    save(strcat(location,'\',sprintf('net%g.mat',i)), 'netTemp');
end

% *** plot training ROC-curve and store training activations ***
if getTrainData
if ~trainNet
    netTemp = networks{i};
end
if plotTrain
    figure
end
for aa=1:4
    Ntrain = NumHSrec.train(aa);
    Ypred = predict(netTemp, getUnbalancedSet(Xtrain{aa},Ntrain));
    YpredTarget = getUnbalancedSet(Ytrain{aa},Ntrain);
    if plotTrain
        subplot(2,2,aa)
    end
    [~,~,~,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        NumHSrec.train(aa),N.segPerPCG,plotTrain);
    if plotTrain
        title(sprintf('AUC training set = %g, pos. %g',AUC.whole,aa))
    end
    CVresults.train.activ{i,aa} = p.murWhole;
end
end

% *** plot validation ROC-curves and store val-set activations ***
if ~trainNet
    netTemp = networks{i};
end
YpredTot = cell(4,1);
YvalTot  = cell(4,1);
if plotVal
    figure
end
for aa=1:4
   % *** get validation data ***
    Ypred = predict(netTemp,Xval{aa});
    if gradedOutput
        YpredTarget = Yval{aa}>=classThr;
    else
        YpredTarget = Yval{aa};
    end
    if plotVal   
        subplot(3,2,aa)
    end
    
    thr = 2;
%     YpredTarget = HSdata.(sprintf('murGrade%g',aa))(CVresults.val.I{i,aa})>=thr;
    [AUC,~,~,~,p,YvalWhole] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        size(Xval{aa},1)/N.segPerPCG,N.segPerPCG,plotVal);
    AUC
    CVresults.val.activ{i,aa} = p.murWhole;
    YvalTot{aa}  = YvalWhole;
    title(sprintf('AUCwhole=%.3g, location=%g',AUC.whole,aa))
    pause(.2)
    AUCmat(i,aa) = AUC.whole;
end

YpredMat = cell2mat(CVresults.val.activ(i,:)');
YvalMat  = cell2mat(YvalTot);
if plotVal
    subplot(3,2,[5,6])
end
AUC = performanceSummaryNeurNet([],YvalMat,YpredMat,...
                                [],N.segPerPCG,true);
title(sprintf('AUC = %g',round(AUC.whole,3)))
pause(.2)
AUCmat(i,5) = AUC.whole;


clearvars Xtrain1234 Ytrain1234 Xval1234 Yval1234

end


%% *** save activations and indeces for CV-sets, and get activation matrix ***
for i=1:Nsplits
    for aa=1:4
        Jtrain = find(CVresults.train.I{i,aa});
        Jval   = find(CVresults.val.I{i,aa});
        CVresults.train.J{i,aa}  = Jtrain;
        CVresults.train.id{i,aa} = ind2id(Jtrain,HSdata);
        CVresults.val.J{i,aa}  = Jval;
        CVresults.val.id{i,aa} = ind2id(Jval,HSdata);
    end
end
% get activation matrix:
CVresults.val.activMat = getZeroPaddedActivMatrix(CVresults.val.activations,...
                                             CVresults.val.J,height(HSdata));

%% ***  Save networks in a cell array ***
networks = cell(8,1);
for i=1:8
    load(sprintf('net%g',i))
    networks{i} = netTemp;
end




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% -------------- TESTING OF OUTPUT ----------------------------------------
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% test 1
net = "regOverTrain";
if net=="reg"
    load 'CVresults_netMurRegAllPos_valStop.mat'
elseif net=="regOverTrain"
    load 'CVresults_netMurRegAllPos_valStop_overTrain.mat'
elseif net=="G2"
    load 'CVresults_netMurG2AllPos_valStop.mat'
elseif net=="G2OverTrain"
    load 'CVresults_netMurG2AllPos_valStop_overTrain.mat'
end

setType = "val";
targetType = "ASgrade";
classThr = 1;
AUCstruc.(net).(setType).(targetType).(sprintf('g%g',classThr)) = zeros(8,5);
targetArray = getTruthCellArray(HSdata,CVresults.(setType).J,targetType,classThr);
close all
plotVal = true;
plotValTot = true;
plotValTotTot = true;
for aa=1:4
    if plotVal
        figure
    end
    for i=1:8
        if targetType=="murmur"
            predVar = sprintf('murGrade%g',aa);
        else
            predVar = targetType;
        end
        Jval  = CVresults.(setType).J{i,aa};
        Ypred = CVresults.(setType).activations{i,aa};
        YpredTarget = targetArray{i,aa};
        if plotVal
            subplot(4,2,i)
        end
        AUC = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                                        [],N.segPerPCG,plotVal);
        AUCstruc.(net).(setType).(targetType).(sprintf('g%g',classThr))(i,aa)...
                        = AUC.whole;
    end
end

for i=1:8
    if targetType=="murmur"
        predVar = sprintf('murGrade%g',aa);
    else
        predVar = targetType;
    end
    Ypred = cell2mat(targetArray(i,:)');
    activ = cell2mat(CVresults.(setType).activations(i,:)');
    AUC = performanceSummaryNeurNet([],Ypred,activ,...
                                    [],N.segPerPCG,plotVal);
    AUCstruc.(net).(setType).(targetType).(sprintf('g%g',classThr))(i,5) = AUC.whole;
end

% *** total AUC for each position ***
if plotValTot
    figure
end
for aa=1:4
    J = cell2mat(CVresults.(setType).J(:,aa));
    if targetType=="murmur"
        predVar = sprintf('murGrade%g',aa);
    else
        predVar = targetType;
    end
    activ = cell2mat(CVresults.(setType).activations(:,aa) );
    target = HSdata.(predVar)(J)>=classThr;
    if plotValTot
        subplot(2,2,aa)
    end
    performanceSummaryNeurNet([],target,activ,...
                              [],N.segPerPCG,plotValTot);
end

% *** total AUC all positions combined ***
if plotValTotTot
    figure
end
J = cell2mat(CVresults.(setType).J(:));
activ = cell2mat(CVresults.(setType).activations(:) );
Ytarget     = cell2mat(targetArray(:));
% target = HSdata.(sprintf('murGrade%g',aa))(J)>=2;
[AUC,X,Y,~,p,YvalWhole] = performanceSummaryNeurNet([],Ytarget,activ,...
                                                    [],N.segPerPCG,plotValTotTot);

%% test 2

close all
load CVresults_netMurRegAllPos_valStop.mat
figure
AUCreg = cell(3,1);
% AUCreg = cell(3,1);
for murThr=1:3
AUCreg{murThr} = zeros(8,4);
for aa=1:4
for i=1:8
    subplot(4,2,i)
    setType = "val";
    Jval  = CVresults.(setType).J{i,aa};
    Ypred = CVresults.(setType).activ{i,aa};
%     YpredTarget = HSdata.ASgrade(Jval)>=1;
    YpredTarget = HSdata.(sprintf('murGrade%g',aa))(Jval)>=murThr;
    [AUC,X,Y,~,p,YvalWhole] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                                                        [],N.segPerPCG,true);
    AUCreg{murThr}(i,aa) = AUC.whole;
end
end
end



