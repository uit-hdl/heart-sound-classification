%% get lengths of each sequence
% Nds = 35;
% Ntrain = height(HSdataTrain);
% lengths = zeros(1,Ntrain);
% for i=1:Ntrain
%     id = HSdataTrain.UNIKT_LOPENR(i);
%     x = downsample(wav2TS(id,1),Nds);
%     lengths(i) = length(x);
% end
% [minLength,i_min] = min(lengths);
% id = HSdataTrain.UNIKT_LOPENR(i_min);
% x_min = downsample(wav2TS(id,1),Nds);
% overLap = floor(segSize*2/3);
% [~,K] = windowSegment(x_min,segSize,overLap)
%% plot some scaleograms
%#ok<*NOPTS>
k=k+1;
id = HSdataTrain.UNIKT_LOPENR(k);
aa = 1;
x0 = wav2TS(id,aa);

% pre process
Nds = 35;
x = downsample(x0,Nds);
segSize   = 3500;
numel(x)
M = getScaleogram(x);
Ncomp = [30,200];
vertSlice = 10:50;
S{1} = imresize(M(vertSlice,1:segSize), Ncomp);
imagesc(M)
imagesc(S{1})
%% generate training data and labes (AR grade are the labels)
fs = 44100;
% murmur threshold:
mg = 2;
% how much to downsample the signal:
Nds = 35;
order = 5;
% auscultation area:
aa = 2;
% segment size to feed the network:
segSize = 3500;
% to which size do we downsize the scaleogram segments:
Ncomp = [30,200];
% number of segments to extract per sample:
Nseg = 3;
% ySlice represents what vertical part of the scaleogram to use:
vertSlice = 10:50;  

% index for cases where there is murmur:
IposTrain = HSdataTrain.ARGRADE_T72>=mg;
JposTrain = find(IposTrain);
IposVal   = HSdataVal.ARGRADE_T72>=mg;
JposVal   = find(IposVal);

Ntrain = height(HSdataTrain);
% overlap of the segments:
overLap = floor(segSize/2);
% Nbalance is number of additional resamplings that is done to balance the
% dataset, so that Nmur==NnoMur:
Nbalance = sum(~IposTrain) - sum(IposTrain);
% size of inflated data set (added samples to account for few positive samples):
NtrainTot = Ntrain + Nbalance;
clearvars Xtrain Ytrain
Xtrain = cell(NtrainTot,Nseg);
Ytrain = zeros(NtrainTot,Nseg);
% Mtrain = cell(NtrainTot,1);

for i=1:NtrainTot
    i 
    if i<=Ntrain
        ii = i;
        id = HSdataTrain.UNIKT_LOPENR(ii);
    else
        % randomly pick a sample of category 1 (done to balance data set)
        ii = randsample(JposTrain,1);
        id = HSdataTrain.UNIKT_LOPENR(ii);
    end

    % get time series:
    x0 = wav2TS(id,aa);
    x  = downsample(x0,Nds);
    
    % get scaleogram and extract overlapping segement indeces:
%     M = getScaleogram(x);
%     Mtrain{i} = M;
    M = Mtrain{i};
    M = M(vertSlice,:);
    [I,K] = windowSegment(x,segSize,overLap);
    % cycle through the segments:
    for k=1:Nseg
        % normalize, downsize, and save segment k:
        Mk = M(:,I{k});
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); % normalize
        MkNorm = imresize(MkNorm, Ncomp);
        Xtrain{i,k} = MkNorm;
        Ytrain(i,k) = IposTrain(ii);
    end
end

Ytrain = categorical(Ytrain);
% reshape to a "vector" form:
Ytrain = reshape(Ytrain,[Nseg*NtrainTot,1]);
Xtrain = reshape(Xtrain,[Nseg*NtrainTot,1]);
%% some zero padding stuff to account for different lengths (not relevant right now...)
% numObservations = numel(Xtrain);
% sequenceLengths = zeros(1,numObservations);
% for i=1:numObservations
%     sequence = Xtrain{i};
%     sequenceLengths(i) = size(sequence,2);
% end
% [sequenceLengths,idx] = sort(sequenceLengths);
% Xtrain = Xtrain(idx);
% Ytrain = Ytrain(idx);
%% TEST SET
Nval = height(HSdataVal);
NsegVal = 2;
Xval = cell(Nval,NsegVal);
Yval = zeros(Nval,NsegVal);
% Mval = cell(Nval,1);

overLap = 0;

for i=1:Nval
    i 
    
    id = HSdataVal.UNIKT_LOPENR(i);
    % get time series:
    x0 = wav2TS(id,aa);
    x  = downsample(x0,Nds);
    
    % get scaleogram and extract overlapping segement indeces:
%     M = getScaleogram(x);
    M = Mval{i};
    M = M(vertSlice,:);
    [I,K] = windowSegment(x, segSize, overLap);
    % cycle through the segments:
    for k=1:NsegVal
        % normalize, downsize, and save segment k:
        Mk = M(:,I{k});
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); % normalize
        MkNorm = imresize(MkNorm, Ncomp);
        Xval{i,k} = MkNorm;
        Yval(i,k) = IposVal(i);
    end 
end

Yval = categorical(Yval);
% reshape to a "vector" form:
Yval = reshape(Yval,[NsegVal*Nval,1]);
Xval = reshape(Xval,[NsegVal*Nval,1]);
%% Define LSTM architecture
numFeatures = height(Xtrain{1}); % 1 timeseries as input per classification
inputSize = numFeatures;

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
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',5,...
    'initialLearnRate',0.001);
%% Train network
net = trainNetwork(Xtrain,Ytrain,layers,options);
%% predict
Ypred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
%%
clear I
mur   = categorical(1);
noMur = categorical(0);
I.pred.mur   = (Ypred==mur)';
I.val.mur    = (Yval==mur)';

acc = sum(Ypred == Yval)./numel(Yval);

NcorrPredAndMur = sum(and(Ypred == Ypred,Ypred==mur));

sens = sum(and(Ypred == Yval,Ypred==mur))/sum(Yval==mur)

spec = sum(and(Ypred==Yval,Ypred==noMur))/sum(Yval==noMur)

%% COMPUTE RUC CURVE
act = activations(net,Xval,'softmax');

I.pred.mur = act(2,:)>.01
xx = 0.00001:.0001:0.99;
xx = (xx);
for i=1:numel(xx)
I.pred.mur = act(2,:)>xx(i)
sens(i) = sum(and(I.pred.mur,I.val.mur))/sum(Yval==mur);
spec(i) = sum(and(~I.pred.mur,~I.val.mur))/sum(Yval==noMur);
end

plot(sens,spec)
ylim([0,1])
ylabel('specificity')
xlabel('sensitivity')
% at best, the classifier achieves a specificity of 92% and a sensitivity
% of 50%. I used 3 segments for each recording, with a 50% overlap for each
% segment, with segment sizes = 3500. The low sensitivity is likely a
% consequence of the fact that there are so few positive cases (murmur
% degree>=3), only 37. Maybe the performance will improve if I use a lower
% murmur threshold.

% for murmur threshold 2, there are 157 cases. The prediction ability is
% about the same as before. I do not think much progress will be made by
% changing the murmur threshold. It is time to try out some segmentation.

% got specificity of 76% and sensitivity of 25% for AR>=2 using the aortic
% position.

% got specificity of 75% and sensitivity of 25% for AR>=2 using the pulonic
% position.