%% get lengths of each sequence
%#ok<*NOPTS>
n = height(HSdataTrain);
Nds = 30;
lengths = zeros(1,n);
for i=1:n
    i
    id = HSdataTrain.UNIKT_LOPENR(i);
    lengths(i) = length(wav2TS(id,1));
end

minLength = floor(min(lengths)/Nds);
%% plot some scaleograms
% k=k+1;
id = HSdataTrain.UNIKT_LOPENR(k);
aa = 1;
x0 = wav2TS(id,aa);

% pre process
Nds = 35;
% passWidth = [25,400];
% x         = preProcessSignal(x0,fs,Nds,passWidth,order);
x = downsample(x0,Nds);
segSize   = 4000;
numel(x)
M = getScaleogram(x);
Ncomp = [30,200];
ySlice = 1:62; 
S{1} = imresize(M(ySlice,1:segSize), Ncomp);
imagesc(M)
imagesc(S{1})
%% generate training data and labes (AR grade are the labels)
n = height(HSdataTrain);
Xtrain = cell(3*n,1);
Ytrain = zeros(3*n,1);

fs = 44100;
% murmur threshold:
mg = 3;
% how much to downsample the signal:
Nds = 35;
order = 5;
% auscultation area:
aa = 1;
% segment size to feed the network:
segSize = 4000;
% to which size do we downsize the scaleogram segments:
Ncomp = [30,200];
% ySlice represents what vertical part of the scaleogram to use:
ySlice = 1:62; 

IposCases = HSdataTrain.Murmur_1_grade_ref_ny_T72>=mg;
JposCases = find(IposCases);
% deltaX represents the additional resampling that is done to balance the
% dataset
deltaX = sum(~IposCases)-sum(IposCases);

n = height(HSdataTrain);
overLap = floor(segSize/2);
nTot = n + deltaX; % size of inflated data set (to account for few positive samples)
Xtrain = cell(nTot,1);
Ytrain = zeros(nTot,1);

for i=1:nTot
    i 
    if i<=n
        id = HSdataTrain.UNIKT_LOPENR(i);
    else
        id = HSdataTrain.UNIKT_LOPENR(randsample(JposCases,1));
    end
        x0        = wav2TS(id,aa);
%         passWidth = [25,400];
%         x         = preProcessSignal(x0,fs,Nds,passWidth,order);
        x = downsample(x0,Nds);
%         segSize   = floor(numel(x)/3);
        
        [I,K] = windowSegment(x,segSize,overLap);
        segs{1} = 1:segSize;
        segs{2} = segSize:2*segSize;
        segs{3} = 2*segSize:3*segSize;
        
        M = getScaleogram(x);
        S{1} = imresize(M(ySlice,1:segSize), Ncomp);
        S{2} = imresize(M(ySlice,segSize:segSize*2), Ncomp);
        
        if 3*segSize>numel(x)
            S{3} = S{2};
        else
            S{3} = imresize(M(ySlice,2*segSize:3*segSize), Ncomp);
        end
        
%         xx{1} = normalise_signal(x(1:segSize));
%         xx{2} = normalise_signal(x(segSize:segSize*2));
%         xx{3} = normalise_signal( x(2*segSize:3*segSize));
    %     x = schmidt_spike_removal(x, fs);
    %     x = normalise_signal(x);
    %     cutoff = 300; % frequenczy in hertz
    %     order = 4;

    %      getScaleogram(x,fs,true);
    %     x = butterworth_low_pass_filter(x,order,cutoff,fs, false);
    %     x = downsample(x,Nds);
        for j = 1:3
            
            Xtrain{3*(i-1) + j,1} = S{j};
            
            if i<=n
                Ytrain(3*(i-1) + j,1) = HSdataTrain.Murmur_1_grade_ref_ny_T72(i)>=mg;
            else
                Ytrain(3*(i-1) + j,1) = (1==1);
            end
            
        end
        
end
Ytrain = categorical(Ytrain);
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
m = height(HSdataVal);
Xval = cell(2*m,1);
Yval = zeros(2*m,1);
for i=1:m
    i
    id = HSdataTrain.UNIKT_LOPENR(i);

    x0        = wav2TS(id,aa);
%     passWidth = [25,400];
%     x         = preProcessSignal(x0,fs,Nds,passWidth,order);
%     segSize   = floor(numel(x)/3);
    x = downsample(x0,Nds);
    
    M    = getScaleogram(x);
    S{1} = imresize(M(ySlice,1:segSize), Ncomp);
    S{2} = imresize(M(ySlice,segSize:2*segSize), Ncomp);
    
    Xval{2*(i-1) + 1,1} = S{1};
    Xval{2*(i-1) + 2,1} = S{2};
    Yval(2*(i-1) + 1,1) = HSdataVal.Murmur_1_grade_ref_ny_T72(i)>=mg;
    Yval(2*(i-1) + 2,1) = HSdataVal.Murmur_1_grade_ref_ny_T72(i)>=mg;
end
Yval = categorical(Yval);

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

%%
act = activations(net,Xval,'softmax');

I.pred.mur = act(2,:)>.028
xx = 0.02:.001:0.4;
for i=1:numel(xx)
sens(i) = sum(and(I.pred.mur,I.val.mur))/sum(Yval==mur);
spec(i) = sum(and(~I.pred.mur,~I.val.mur))/sum(Yval==noMur);
end

plot(sens)
