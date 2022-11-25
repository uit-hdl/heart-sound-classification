load 


%% Train network, or get training or validation set
%#ok<*NOPTS>
load nodesAll
nodes = nodesAll;

% storage variables:
activations = cell(1,4);
JtrainCell = cell(1,4);
JvalCell   = cell(1,4);
Xtrain = cell(1,4);
Ytrain = cell(1,4);
Xval = cell(1,4);
Yval = cell(1,4);

% ¤¤ CHOOSE IF TRAIN NETWORK ¤¤
trainNet = false;
% ¤¤ CHOOSE IF GET TRAINING DATA ¤¤
getTrainingData = true;
% ¤¤ CHOOSE WHETHER OR NOT TO USE VALIDATION ACCURACY DURING TRAINING ¤¤
checkValAccuracy = false;
for aa=1:4
    aa 
    % define mother dataframe:
    data0    = HSdata;
    % remove rows corresponding to noissy observations:
    Jclean   = find(data0.(sprintf('MURMUR_%gNOISE_REF_T72',aa))==0);
    data0    = data0(Jclean,:);
    % get segmentation info for 
    nodes0.loc      = nodes.loc(Jclean,:);
    nodes0.state    = nodes.state(Jclean,:);
    nodes0.seglines = nodes.seglines(Jclean,:);
    Jtrain = getNewAdress(HSdata,union(Jtrain0,Jval0),Jclean);
    Jval   = getNewAdress(HSdata,Jtest0,Jclean);
    
    % define the labels
    murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    Y0    = data0.(murSt);
    N0    = height(data0);

    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
    
    posThr = 1;
    if trainNet
        balanceTrain = true;
    else
        balanceTrain = false;
    end
    balanceVal   = false;
    Ncomp        = [30,200];
    
    if getTrainingData
       [Xtrain{aa},Ytrain{aa}] = genTrainOrValSet(data0,Y0,Jtrain,N,aa,...
                                        balanceTrain,nodes0,Ncomp,[],posThr);
    end
   [Xval{aa},Yval{aa}] = genTrainOrValSet(data0,Y0,Jval,N,aa,...
                                    balanceVal,nodes0,Ncomp);
                                
    JtrainCell{1,aa} = Jclean(Jtrain);
    JvalCell{1,aa}   = Jclean(Jval);
                                
end

% combine into 1 training set and validation set:
if trainNet
    Xtrain1234 = [Xtrain{1};Xtrain{2};Xtrain{3};Xtrain{4}];
    Ytrain1234 = [Ytrain{1};Ytrain{2};Ytrain{3};Ytrain{4}];
end
Xval1234 = [Xval{1};Xval{2};Xval{3};Xval{4}];
Yval1234 = [Yval{1};Yval{2};Yval{3};Yval{4}];

%%% train %%%
if trainNet
    % *** define architecture ***
    numFeatures = height(Xtrain{1}{1}); % same as input size...
    inputSize   = numFeatures; % number of timeseries to take as input
    numHiddenUnits = 50;
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
            lstmLayer(numHiddenUnits)
            lstmLayer(numHiddenUnits,'OutputMode','last')
            dropoutLayer(.5)
            fullyConnectedLayer(30)
            reluLayer
            fullyConnectedLayer(1)
            regressionLayer];
    end
    
    % *** specify training options ***
    % ¤¤ CHOOSE MAX NUMBER OF EPOCHS ¤¤
    maxEpochs = 50;
    % ¤¤ SET INITIAL LEARN RATE AND LEARN-RATE-DROP-PERIOD ¤¤
    initLearnRate = 0.002;
    LearnRateDropPeriod = 5;
    miniBatchSize = 2^5;
    validationFrequency = 500;%floor(numel(Ytrain)/miniBatchSize);
    if checkValAccuracy
        validationData = {XvalBal1234,YvalBal1234};
    else
        validationData = [];
    end
    options = trainingOptions('adam', ...
        'ExecutionEnvironment','cpu', ...
        'MaxEpochs',maxEpochs, ...
        'MiniBatchSize',miniBatchSize, ...
        'GradientThreshold',1, ...
        'Plots','training-progress',...
        'LearnRateSchedule','piecewise', ...
        'LearnRateDropFactor',0.5, ...
        'LearnRateDropPeriod',5, ...
        'initialLearnRate',initLearnRate, ...
        'ValidationData',validationData, ...
        'ValidationFrequency',validationFrequency, ...
        'LearnRateDropPeriod',LearnRateDropPeriod, ...
        'initialLearnRate',initLearnRate, ...
        'OutputFcn',@(info)myCostumOutputFcn(info,270*60,inf));
    
    % *** train and save network ***
    net = trainNetwork(Xtrain1234,Ytrain1234,layers,options);
end

%% training ROC-curve
figure
for aa=1:4
    Ypred = predict(net,Xtrain{aa},'MiniBatchSize',miniBatchSize);
    % ¤¤ CHOOSE PREDICITON TARGET:
    YpredTarget = Ytrain{aa}>=2;
    subplot(2,2,aa)
    [AUC,~,~,~,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        size(Xtrain{aa},1)/N.segPerPCG,N.segPerPCG,true,1);
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
    activations{aa} = p.murWhole;
    pause(.6)
end
%% test set predictions
% clearvars Xtrain1234 Ytrain1234 Xval1234 Yval1234

YpredCell = cell(4,1);
YvalCell  = cell(4,1);
figure
for aa=1:4
    Ypred  = predict(net,Xval{aa},'MiniBatchSize',miniBatchSize);
    % ¤¤ CHOOSE PREDICITON TARGET:
%     YpredTarget = Yval{aa}>=2;
    YpredTarget = HSdata.ASGRADE_T72(JvalCell{1,aa})>=1;
    subplot(2,2,aa)
    [AUC,~,~,~,p,YvalWhole] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        size(Xval{aa},1)/N.segPerPCG,N.segPerPCG,true);
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
    pause(.1)
    YpredCell{aa} = p.murWhole;
    YvalCell{aa}  = YvalWhole;
end
Ypred    = cell2mat(YpredCell);
YvalAll  = cell2mat(YvalCell);

YpredTarget = YvalAll;
figure
[AUC] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                                  [],N.segPerPCG,true);
                              
%% simpler way of doing it:

J = cell2mat(CVresults.val.J')

classThr = 2;
Yval = cell(1,4);
for aa=1:4
    J = CVresults.val.J{1,aa}
    targetVar = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    Yval{1,aa} = HSdata.(targetVar)(J)>=classThr;
end

YpredTarget = cell2mat(Yval');
Ypred = cell2mat(CVresults.val.activations');
[AUC] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                                  [],[],true);
legend('murmur-cutoff=1 (AUC=0.936)', 'murmur-cutoff=2 (AUC=0.966)')


