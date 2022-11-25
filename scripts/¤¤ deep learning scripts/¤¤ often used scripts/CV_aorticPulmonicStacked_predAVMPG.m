% This script trains networks (any target) on the
% cross-validation-partitions, using data from all four positions.

% Load data on CV-partitioning and signal segmentation:
load nodesAll.mat
load CVpartitions.mat
nodes0 = nodesAll;

Nsplits = 8;

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
targetType = "avmeanpg";
% ¤¤ CHOOSE IF RESAMPLE TO BALANCE DATA ¤¤
balanceTrain = true;
balanceVal   = true;
% ¤¤ CHOOSE IF TRAINING OUTPUT VARIABLE IS GRADED OR BINARY ¤¤
gradedOutput = true;
binaryOutput = ~gradedOutput;

if not(gradedOutput)
    % ¤¤ IF BINARY OUTPUT: CHOOSE CLASS SPERATION THRESHOLD ¤¤
    classThr = 2;
    posThr = classThr;
elseif gradedOutput
    % ¤¤ IF GRADED OUTPUT: SET THRESHOLD FOR RESAMPLING CLASS ¤¤
    posThr = 6.0;
    classThr = 10.0; % class threshold for plotting
end
% --- note: resampling is done for all Y >= posThr ---
% ¤¤ PLOT SETTINGS ¤¤
plotTrain = false;
plotVal = true;
close all


Info = table(trainNet,getTrainData,balanceTrain,balanceVal,binaryOutput,...
             loadNetsName)

P = [1,2];
Npos = numel(P);

% define some result storage variables:
CVresults.val.I     = cell(Nsplits,Npos);
CVresults.val.activ = cell(Nsplits,Npos);
CVresults.train.I     = cell(Nsplits,Npos);
CVresults.train.activ = cell(Nsplits,Npos);
CVresults.valTot.I   = cell(Nsplits,1);
CVresults.trainTot.I = cell(Nsplits,1);

% keep track of number of recordings for each validation set and position:
NumHSrec.val   = zeros(Npos,1);
NumHSrec.train = zeros(Npos,1);
AUCmat = zeros(Nsplits,Npos); % the fifth column is for the AUC corresponding to all positions

for i=1:8
disp(i) %#ok<*NOPTS>
Xtrain = cell(Npos,1);
Ytrain = cell(Npos,1);
XvalBal = cell(Npos,1);
YvalBal = cell(Npos,1);
Xval   = cell(Npos,1);
Yval   = cell(Npos,1);

% get index for rows corresponding to the training and validation sets:
ItrainRows = CVpartitions.train.I{i};
IvalRows   = CVpartitions.val.I{i};
IavmpgHasValue = isval(HSdata.avmeanpg);
% $$$ $$$ $$$-- get validation and/or training data --$$$ $$$ $$$
Iclean = or(HSdata.noise1==0,HSdata.noise2==0);
Iclean = and(Iclean,IavmpgHasValue);

for aa=1:2
    disp(aa)
    % get index of clean recordings for position aa:
    % get index for training and validation set observations of position aa:
    Jnoise = find(HSdata.(sprintf('noise%g',aa))==1);
    Itrain = and(ItrainRows,Iclean);
    Ival   = and(IvalRows,Iclean);
    Jtrain = find(Itrain);
    Jval   = find(Ival);
    
    targetStr = targetType;


    if gradedOutput
        Y0 = HSdata.(targetStr);
    elseif binaryOutput
        Y0 = HSdata.(targetStr)>=classThr;
    end
    
    % segment extraction settings:
    N.ds = 20; % how much to downsample signal
    N.cs = 4;  % How many heartbeats per segment
    N.os = 2;  % by how many heartbeats do the segments overlap
    N.segPerPCG = 6; % how many segments to extract per recording
    Ncomp = [13,200]; % MFCC of segments are reshaped to these dimensions
    
   if aa==1
       JresampTrain = [];
       JresampVal = [];
   end
    
    if trainNet || getTrainData
       % *** get training data ***
       [Xtrain{aa},Ytrain{aa},~,~,JresampTrain] = genTrainOrValSet_new(HSdata,Y0,{Jnoise,Jtrain},N,aa,...
                                        balanceTrain,nodes0,Ncomp,[],posThr,JresampTrain);
    end
    
   % *** get validation data ***
   [XvalBal{aa},YvalBal{aa},Xval{aa},Yval{aa},JresampVal] = genTrainOrValSet_new(HSdata,Y0,{Jnoise,Jval},N,aa,...
                                    balanceVal,nodes0,Ncomp,[],posThr,JresampVal);
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

% Stack data from all pos. into one training and validation set to feed
% to training algorithm:
NtrainSamples = numel(Xtrain{1});
XtrainStacked = cell(NtrainSamples,1);
for k=1:NtrainSamples
    XtrainStacked{k} = [Xtrain{1}{k};Xtrain{2}{k}];
end
YtrainStacked = Ytrain{1};

NvalBalSamples = numel(XvalBal{1});
XvalBalStacked = cell(NvalBalSamples,1);
NvalSamples = numel(Xval{1});
XvalStacked = cell(NvalSamples,1);
for k=1:NvalBalSamples
    XvalBalStacked{k} = [XvalBal{1}{k};XvalBal{2}{k}]; 
end
for k=1:NvalSamples
    XvalStacked{k} = [Xval{1}{k};Xval{2}{k}]; 
end
YvalBalStacked = YvalBal{1};
YvalStacked = Yval{1};

% $$$ $$$ $$$-- train network --$$$ $$$ $$$ 
if trainNet
    % *** define network architecture ***
    numFeatures = height(XtrainStacked{1}); % same as input size...
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
    checkValAccuracy = true;
    miniBatchSize = 2^5;
    validationFrequency = floor(numel(XtrainStacked)/miniBatchSize);
    if checkValAccuracy
        validationData = {XvalBalStacked,YvalBalStacked};
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
        'initialLearnRate',0.0005, ...
        'ValidationData',validationData, ...
        'ValidationFrequency',validationFrequency, ...
        'LearnRateDropPeriod',5, ...
        'OutputFcn',@(info)myCostumOutputFcn(info));
    
    % *** train and save network ***
    % ¤¤ CHOOSE IF INITIATE WITH WEIGHTS FROM PRETRAINED NET ¤¤
    initTrainingWithNets = false;
    initNetsName = nan;
    if initTrainingWithNets
        initNetsName = "networksCVnoNoiseDrOutMurRegAllPosValStop.mat";
        initNets = load(initNetsName);
        initNet = initNets.networks{i};
    end
    InfoTrain = table(checkValAccuracy,initTrainingWithNets,initNetsName)
    
    if initTrainingWithNets
        netTemp = trainNetwork(XtrainStacked,YtrainStacked,initNet.Layers,options);
    else
        netTemp = trainNetwork(XtrainStacked,YtrainStacked,layers,options);
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

% *** get validation data ***
Ypred = predict(netTemp,XvalStacked);
if gradedOutput
    YpredTarget = YvalStacked>=classThr;
else
    YpredTarget = YvalStacked;
end

YpredTarget = HSdata.avmeanpg(CVresults.valTot.I{i})>=10;
[AUC,~,~,~,p,YvalWhole] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
    size(XvalStacked,1)/N.segPerPCG,N.segPerPCG,plotVal);
AUC
CVresults.val.activ{i,aa} = p.murWhole;
YvalTot{aa}  = YvalWhole;
title(sprintf('AUCwhole=%.3g, location=%g',AUC.whole,aa))
pause(.2)
AUCmat(i,aa) = AUC.whole;

YpredMat = cell2mat(CVresults.val.activ(i,:)');
YvalMat  = cell2mat(YvalTot);
if plotVal
    subplot(3,2,[5,6])
end
AUC = performanceSummaryNeurNet([],YvalMat,YpredMat,...
                                [],N.segPerPCG,true);
pause(.2)
AUCmat(i,5) = AUC.whole;

clearvars Xtrain1234 Ytrain1234 Xval1234 Yval1234

end
AUCcont{end+1} = AUCmat


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
