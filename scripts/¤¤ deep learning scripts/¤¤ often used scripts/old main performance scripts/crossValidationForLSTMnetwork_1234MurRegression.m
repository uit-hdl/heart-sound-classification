% In this script, I train
load nodes

Nsplits = 8;

% ¤¤ CHOOSE WHETHER TO TRAIN OR LOAD PRETRAINED NETWORK FOR TESTING ¤¤
trainNet = true;
if trainNet
    networks = cell(1,8);
else
    load networksCVnoNoiseDrOutRegaa1234.mat
end

activations = cell(Nsplits,4);
AUCmat      = zeros(Nsplits,4);

for i=1:Nsplits
i %#ok<*NOPTS>
Xtrain = cell(4,1);
Ytrain = cell(4,1);
XvalBal = cell(4,1);
YvalBal = cell(4,1);
Xval   = cell(4,1);
Yval   = cell(4,1);

for aa=1:4
    % get clean data:
    data0    = HSdataTrain; % Use training set data for CV
    % get 
    Jnonoise = data0.(sprintf('MURMUR_%gNOISE_REF_T72',aa))==0;
    data0    = data0(Jnonoise,:);
    nodes0.loc      = nodes.loc(Jnonoise,:);
    nodes0.state    = nodes.state(Jnonoise,:);
    nodes0.seglines = nodes.seglines(Jnonoise,:);

    Y0   = data0.(murSt)>=2;
    N0   = height(data0);
    Nval = floor(N0/Nsplits);
    Jall = 1:N0;

    Jval   = (i-1)*Nval+1:i*Nval;
    Jtrain = setdiff(Jall,Jval);

    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
    
    balanceTrain = true;
    balanceVal   = true;
    Ncomp = [30,200];

   [Xtrain{aa},Ytrain{aa}] = genTrainOrValSet(data0,Y0,Jtrain,N,aa,...
                                    balanceTrain,nodes0,Ncomp);
   [XvalBal{aa},YvalBal{aa},Xval{aa},Yval{aa}] = genTrainOrValSet(data0,Y0,Jval,N,aa,...
                                    balanceVal,nodes0,Ncomp);
end
% combine into 1 training set and validation set:
Xtrain1234 = [Xtrain{1};Xtrain{2};Xtrain{3};Xtrain{4}];
Ytrain1234 = [Ytrain{1};Ytrain{2};Ytrain{3};Ytrain{4}];
XvalBal1234   = [XvalBal{1};XvalBal{2};XvalBal{3};XvalBal{4}];
YvalBal1234   = [YvalBal{1};YvalBal{2};YvalBal{3};YvalBal{4}];

%%% train %%%
if trainNet
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
    % specify training options
    % ¤¤ CHOOSE MAX NUMBER OF EPOCHS:
    maxEpochs = 100;
    % ¤¤ SET INITIAL LEARN RATE AND LEARN-RATE-DROP-PERIOD ¤¤
    initLearnRate = 0.002;
    LearnRateDropPeriod = 5;
    miniBatchSize = 2^5;
    validationFrequency = floor(numel(Ytrain)/miniBatchSize);
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
    
    % train and save network:
    networks{i} = trainNetwork(Xtrain1234,Ytrain1234,layers,options);
end

Nloc = 1;
figure
for aa=1:4
    Ypred = predict(networks{i},Xval{aa},'MiniBatchSize',miniBatchSize);
    activations{i,aa} = Ypred;
    % ¤¤ CHOOSE PREDICITON TARGET:
    YpredTarget = Yval{aa};
    subplot(2,2,aa)
    [AUC,X,Y,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        size(Xval{aa},1)/N.segPerPCG,N.segPerPCG,true,Nloc);
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
    pause(.1)
    if aa==1
    AUCmat(end+1,aa) = AUC.whole; %#ok<*SAGROW>
    else
    AUCmat(end,aa) = AUC.whole; 
    end
end
clearvars Xtrain1234 Ytrain1234 Xval1234 Yval1234

end



