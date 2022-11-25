dataModel = HSdata;
Ndata = height(dataModel);
Ntrain = floor(Ndata*3/4);
Nval = Ndata-Ntrain;

Jall = 1:Ndata;
Jtrain = randsample(Jall,Ntrain);
Jval = setdiff(Jall,Jtrain);

Y0 = dataModel.ARGRADE_T72>=3;

Jpos = find(Y0>0);
Npos = sum(Y0>0);
Nneg = Ndata - Npos;
NposResample = Nneg - Npos;
JposResample = randsample(Jtrain,NposResample);

Ndiff = Npos;
X0 = [dataModel.Murmur_1_grade_ref_ny_T72,dataModel.Murmur_3_grade_ref_ny_T72...
     dataModel.Murmur_3_grade_ref_ny_T72,dataModel.Murmur_4_grade_ref_ny_T72];

Xtrain = X0(Jtrain,:);
Ytrain = Y0(Jtrain,:);

Xval = X0(Jval,:);
Yval = Y0(Jval,:);
Ytrain = categorical(Ytrain);

numFeatures = width(Xtrain);
numClasses = 2;

layers = [ ...
    featureInputLayer(numFeatures)
    fullyConnectedLayer(20)
    reluLayer
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

maxEpochs = 70;
miniBatchSize = 2^4;
options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Plots','training-progress',...
    'LearnRateSchedule','piecewise',...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',5, ...
    'initialLearnRate',0.002);

net = trainNetwork(Xtrain,Ytrain,layers,options);

Ypred = predict(net,Xval,'MiniBatchSize',miniBatchSize);
