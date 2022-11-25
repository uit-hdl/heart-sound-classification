murAmp = .2;
L = 12000; % length of a full sample
T = [3000,3000]; %  periods for class 1 and 2
loc = [.1 .4];
offsets = randi(L,2,1);
x1 = simulHSwithVarFreq(L,T(1),murAmp,loc,offsets(1));
x2 = simulHSwithVarFreq(L,T(2),0,loc,offsets(2));

Ncomp = [15,110];
ySlice = 15:68;
clf
subplot(211)
plot(x1)
scg = getScaleogram(x1,1,true);
scg = imgaussfilt(scg,1);
scgCompressed = imresize(scg(ySlice,:), Ncomp);
imagesc((scgCompressed))
subplot(212)
plot(x2)
scg = getScaleogram(x2,1,true);
scg = imgaussfilt(scg,1);
scgCompressed = imresize(scg(ySlice,:), Ncomp);
imagesc((scgCompressed))
%% SIMULATE TRAINING DATA
Ntrain = 1400;
L = 12000; % length of a full "recording"
T = [3000, 3000]; %  periods for class 1 and 2
loc = [.1 .4];
overLap = 0;
Nseg = 1;   
Xtrain = cell(Ntrain,Nseg);
Ytrain = zeros(Ntrain,Nseg);
% segSize = floor(L/Nseg);
segSize = T(1)*4;
murAmp1 = .2;
murAmp2 = 0;

for i=1:Ntrain
    i
    if i<Ntrain/2
        % simulate
        r = randi(L,1);
        x = simulHSwithVarFreq(L,T(1),murAmp1,loc,r);
        
        
        M = getScaleogram(x);
        M = imresize(M(ySlice,:), Ncomp);
        [I,K] = windowSegment(x,segSize,overLap);
        for k=1:Nseg
            % normalize scaleogram segment:
            xx = (M-mean(M,'all'))/std(M,0,'all');
            scgCompressed = imresize(xx, Ncomp);
            Xtrain{i,k} = xx;
            Ytrain(i,k) = 1;
        end
        
    else
        murAmp = 0;
        r = randi(L,1);
        x = simulHSwithVarFreq(L,T(2),murAmp2,loc,r);
        M = getScaleogram(x);
        M = imresize(M(ySlice,:), Ncomp);
        [I,K] = windowSegment(x,segSize,0);
        
        for k=1:Nseg
            xx = (M-mean(M,'all'))/std(M,0,'all'); % normalize
            scgCompressed = imresize(xx, Ncomp);
            Xtrain{i,k} = xx;
            Ytrain(i,k) = 2;
        end
    end
end

Ytrain = categorical(Ytrain);
% reshape
Ytrain = reshape(Ytrain,[Nseg*Ntrain,1]);
Xtrain = reshape(Xtrain,[Nseg*Ntrain,1]);
%% plot to check that all went well
clf
a = 21;
b = floor(length(Ytrain)*0.1/5);
subplot(221)
imagesc(Xtrain{a});
subplot(222)
imagesc(Xtrain{a+1});
subplot(223)
imagesc(Xtrain{b});
subplot(224)
imagesc(Xtrain{b+1});

%% SIMULATE VALIDATION DATA
Nval = 200;
Xval = cell(Nval,Nseg);
Yval = zeros(Nval,Nseg);
for i=1:Nval
    i %#ok<*NOPTS>
    if i<Nval/2
        r = randi(L,1);
        x = simulHSwithVarFreq(L,T(1),murAmp1,loc,r);
        M = getScaleogram(x);
        M = imresize(M(ySlice,:), Ncomp);
        [I,K] = windowSegment(x,segSize,overLap);
        for k=1:Nseg
            xx = (M-mean(M,'all'))/std(M,0,'all'); % normalize
            scgCompressed = imresize(xx, Ncomp);
            Xval{i,k} = xx;
            Yval(i,k) = 1;
        end
        
    else
        murAmp = 0;
        r = randi(L,1);
        x = simulHSwithVarFreq(L,T(2),murAmp2,loc,r);
        M = getScaleogram(x);
        M = imresize(M(ySlice,:), Ncomp);
        [I,K] = windowSegment(x,segSize,0);
        
        for k=1:Nseg
            xx = (M-mean(M,'all'))/std(M,0,'all'); % normalize
            scgCompressed = imresize(xx, Ncomp);
            Xval{i,k} = xx;
            Yval(i,k) = 2;
        end
    end
end

Yval = categorical(Yval);
Yval = reshape(Yval,[Nseg*Nval,1]);
Xval = reshape(Xval,[Nseg*Nval,1]);

clf
subplot(211)
imagesc(Xval{2})
subplot(212)
imagesc(Xval{end-1})
%% Define LSTM architecture
numFeatures = size(Xtrain{1},1); % 10 timeseries as input per classification
inputSize = Ncomp(1);

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
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',5,...
    'initialLearnRate',0.001);
%% Train network
net = trainNetwork(Xtrain,Ytrain,layers,options);
%% predict
YPred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
YPredTrain = classify(net,Xtrain,'MiniBatchSize',miniBatchSize);
%%
accTrain = sum(YPredTrain == Ytrain)./numel(Ytrain)
acc = sum(YPred == Yval)./numel(Yval)
(sum(mode(reshape(YPred,[Nval,Nseg]),2)==mode(reshape(Yval,[Nval,Nseg]),2)))/Nval
% I downsized the scaleograms to size [10,60] (a huge dimension reduction)
% which resulted in very pixuly images. Then I used considered each row of
% this pixulated, reduced form as a feature, so I got 10 feature for each
% sample, and each feature was a time series of length 60. The LSTM
% classifier reached 100% training accuracy in a matter of seconds, and
% reached 100% accuracy on the training set. Clearly, this operation made
% an enormeous difference, as the classifier was completely unable to train
% when I fed it segments (corresponding to cycles) of data as input. 

% reducing the murmur amplitude to 0.4 resulted in slower training time (30
% seconds), but 100% accuray on the test set. The murmur is now getting
% somewhat hard to see.

% Murmur amplitude is now 0.2, and getting quite hard to see. The
% classifier reaches 93% test accuracy, but is getting significantly slower
% to learn. Increasing the learn rate drop period helped it not get stuck.

% follow up from above: I INCREASED the RESOLUTION of the SCALEOGRAM up to
% [15,70]. This did not alter the results significantly. the test error is
% 0.94 and the training error is now 0.996, so we are beginning to see some
% signs of overtraining. It appears that there was MORE OVERTRAINING with
% HIGHER RESOLUTION. This makes sense, as it is harder to overtrain with
% lower dimensional input, as the individual samples loose much of their
% identity due to the averaging operation that takes place when the images
% is downsized.

% follow up from above: running two segments instead of 1 does not seem to
% hurt performance. Running all 4 segments at once does not hurt
% performance.

% When the data has random initial phase (murAmp=0.2), and I do not account
% for this (I feed the whole reduced image, without segmentation) the
% algorithm performance drops significantly from ~93-96% test accuracy to
% ~86% test accuracy. It does manage to perform decently and converge
% withing reasonable time, but the task has clearly become much more
% difficult. Is it possible to remedy this by increasing the amount of
% data? Try doubling the amount of data to see what happens...

% follow up to above: Doubling the amount of training data dramatically
% increased performance, from 85% -> 99% test accuracy. It appears that
% that more data can compensate for lack of segmentation.

% act = activations(net,Xval,'softmax')

