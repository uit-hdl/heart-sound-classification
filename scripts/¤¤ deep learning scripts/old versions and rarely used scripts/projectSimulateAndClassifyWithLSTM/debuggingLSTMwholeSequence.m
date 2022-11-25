Ntrain = 3000;
Xtrain = cell(Ntrain,1);
Ytrain = zeros(Ntrain,1);
L = 1200;
T = [300,400];
loc = [.1 .3 .45];
for i=1:Ntrain
    
    if i<Ntrain/2
        murAmp = 40;
        Xtrain{i} = normalize(simulHeartSound(murAmp,T(1),L,loc));
        Ytrain(i) = 1;
    else
        murAmp = 0;
        Xtrain{i} = normalize(simulHeartSound(murAmp,T(2),L,loc));
        Ytrain(i) = 2;
    end
end

Ytrain = categorical(Ytrain);
clf
subplot(211)
plot(Xtrain{1})
subplot(212)
plot(Xtrain{2000})
%% test set
Nval  = 500;
Xval = cell(Nval,1);
Yval = zeros(Nval,1);
for i=1:Nval
    
    if i<Nval/2
        murAmp = 40;
        Xval{i} = normalize(simulHeartSound(murAmp,T(1),L,loc));
        Yval(i) = 1;
    else
        murAmp = 0;
        Xval{i} = normalize(simulHeartSound(murAmp,T(2),L,loc));
        Yval(i) = 2;
    end
end

Yval = categorical(Yval);
clf
subplot(211)
plot(Xval{1})
subplot(212)
plot(Xval{300})

%
% numObservations = numel(Xtrain);
% sequenceLengths = zeros(1,numObservations);
% for i=1:numObservations
%     sequence = Xtrain{i};
%     sequenceLengths(i) = size(sequence,2);
% end
% [sequenceLengths,idx] = sort(sequenceLengths);
% Xtrain = Xtrain(idx);
% Ytrain = Ytrain(idx);
%% Define LSTM architecture
numFeatures = 1; % 1 timeseries as input per classification
inputSize = 1;

numHiddenUnits = 100;
numClasses = 2;

layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];

%% specify training options
maxEpochs = 70;
miniBatchSize = 2^5;

options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Verbose',false, ...
    'Plots','training-progress',...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.2, ...
    'LearnRateDropPeriod',1,...
    'initialLearnRate',0.001);
%% Train network
net = trainNetwork(Xtrain,Ytrain,layers,options);
%% predict
YPred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
%%
acc = sum(YPred == Yval)./numel(Yval)
    