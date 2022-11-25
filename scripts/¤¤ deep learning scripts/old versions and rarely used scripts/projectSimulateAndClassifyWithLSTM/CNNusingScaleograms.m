murAmp = 10;
L = 2000; % length of a full sample
T = [300,300]; %  periods for class 1 and 2
loc = [.1 .3 .45];
x1 = simulHeartSound(10,T(2),L,loc);
x2 = simulHeartSound(0,T(2),L,loc);
clf
subplot(211)
plot(x1)
subplot(212)
plot(x2)
%%
Ntrain = 2000;
L = 2000; % length of a full sample
T = [300,300]; %  periods for class 1 and 2
loc = [.1 .3 .45];
overLap = 0;
Nseg = 4;
Xtrain = cell(Ntrain,Nseg);
Ytrain = zeros(Ntrain,Nseg);
% segSize = floor(L/Nseg);
segSize = T(1);
murAmp1 = 20;
murAmp2 = 0;

for i=1:Ntrain
    
    if i<Ntrain/2
        x = simulHeartSound(murAmp1,T(1),L,loc);
        [I,K] = windowSegment(x,segSize,overLap);
        for k=1:Nseg
            Xtrain{i,k} = normalise_signal(x(I{k}));
            Ytrain(i,k) = 1;
        end
        
    else
        murAmp = 0;
        x = simulHeartSound(murAmp2,T(2),L,loc);
        [I,K] = windowSegment(x,segSize,0);
        
        for k=1:Nseg
            Xtrain{i,k} = normalise_signal(x(I{k}));
            Ytrain(i,k) = 2;
        end
    end
end

Ytrain = categorical(Ytrain);
Ytrain = reshape(Ytrain,[Nseg*Ntrain,1]);
Xtrain = reshape(Xtrain,[Nseg*Ntrain,1]);
%%
clf
a = 21;
b = 3121;
subplot(221)
getScaleogram(Xtrain{a},1,true);
subplot(222)
getScaleogram(Xtrain{a+1},1,true);
subplot(223)
getScaleogram(Xtrain{b},1,true);
subplot(224)
getScaleogram(Xtrain{b+1},1,true);

%% test set
Nval = 200;
Xval = cell(Nval,Nseg);
Yval = zeros(Nval,Nseg);

for i=1:Nval
    
    if i<Nval/2 
        x = simulHeartSound(murAmp1,T(1),L,loc);
        [I,K] = windowSegment(x,segSize,overLap);
        for k=1:Nseg
            Xval{i,k} = normalise_signal(x(I{k}));
            Yval(i,k) = 1;
        end
        
    else
        x = simulHeartSound(murAmp2,T(2),L,loc);
        [I,K] = windowSegment(x,segSize,0);
        
        for k=1:Nseg
            Xval{i,k} = normalise_signal(x(I{k}));
            Yval(i,k) = 2;
        end
    end
end

Yval = categorical(Yval);
Yval = reshape(Yval,[Nseg*Nval,1]);
Xval = reshape(Xval,[Nseg*Nval,1]);

close all
figure()
subplot(211)
plot(Xval{2})
subplot(212)
plot(Xval{500})

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
    imageInputLayer([28 28 1])
    convolution2dLayer(5,20)
    reluLayer
    maxPooling2dLayer(2,'Stride',2)
    fullyConnectedLayer(10)
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
    'LearnRateDropFactor',0.4, ...
    'LearnRateDropPeriod',2,...
    'initialLearnRate',0.001);
%% Train network
net = trainNetwork(imdsTrain,layers,options);
%% predict
YPred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
%%
acc = sum(YPred == Yval)./numel(Yval)