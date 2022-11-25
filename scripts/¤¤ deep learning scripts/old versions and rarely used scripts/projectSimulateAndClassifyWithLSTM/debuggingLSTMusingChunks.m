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
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(30)
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
    'LearnRateDropFactor',0.4, ...
    'LearnRateDropPeriod',2,...
    'initialLearnRate',0.001);
%% Train network
net = trainNetwork(Xtrain,Ytrain,layers,options);
%% predict
YPred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
%%
acc = sum(YPred == Yval)./numel(Yval)
% conclusion: when the segments become too small, the performance of the
% network starts decreasing. There is clearly a sweet spot in terms of
% segment length.

% conclusion 2: segmenting the data made a huge difference. After removing
% the difference in period lenght, the method of using segments of
% arbitrary lengths completely failed to produce a result beyond random
% luck. When using the true cycle length T to define the segment length
% however, the performance almost immediately went up to 97%, and the
% training was fast. Even when reducing the murmur amplitude 20, the
% training managed to produce good results (93% accuracy). Conclusion:
% Segmentation MATTERS. Note: when the murmur amplitude was set to 20, it
% took the optimization algorithm 4 minutes before it finally latched onto
% a fueature that it could use. When it did latch, it only took a few
% seconds before it settled around ~90% accuracy. It appears that the more
% subtle the feature is, the longer the algorithm is stuck searching for
% the feature that makes the difference.

