% find the number of segments available for each sample

n_samples = zeros(1,height(HSdataTrain));
for i=1:height(HSdataTrain)
    n_samples(i) = numel(segLines{i,1});
end
% find the smallest number of cycles in a sample:
[ncycle_min, imin] = min(n_samples)
ncycle_min = ncycle_min-1

xminCycles = wav2TS(HSdataTrain.UNIKT_LOPENR(imin),aa);
xminCycles = downsample(xminCycles,20)
getScaleogram(xminCycles,imin,true);
% each sample contains a minimum of 4 cycles
%% plot some scaleograms
%#ok<*NOPTS>
load 'segLines.mat'
% k  = k+1;
k=1
id = HSdataTrain.UNIKT_LOPENR(k);
aa = 1;
x0 = wav2TS(id,aa);

% pre process
Nds = 20;
x = downsample(x0,Nds);
NcyclesPerSeg = 5;
clf
M = getScaleogram(x,1,false);
Ncomp = [30,200];
vertSlice = 15:58;
Nos = 2;
Ncs = 4;
NsegDesired = 12;
L = getSegments(segLines{6,1},Nos,Ncs,NsegDesired);
ind_seg = 1;
S{1} = imresize(M(vertSlice,L(ind_seg,1):L(ind_seg,2)), Ncomp );
imagesc(M)
imagesc(S{1})
%% generate training data and labels
load 'segLines.mat'
fs = 44100;
% murmur threshold:
mg = 2;
% how much to downsample the signal:
Nds = 20;
% auscultation area:
aa = 1;
% number of cycles per segment:
Ncs = 4;
% cycle overlap for each consecutive pair of segments:
Nos = 2;
% number of segments to extract from each sample:
NsegDesired = 3;
% to which size do we downsize the scaleogram segments:
Ncomp = [30,200];
% ySlice represents what vertical part of the scaleogram to use:
vertSlice = 10:53;  

% index for cases where there is murmur:
IposCases = HSdataTrain.Murmur_1_grade_ref_ny_T72>=mg;
JposCases = find(IposCases);

Ntrain = height(HSdataTrain);
% Nbalance is number of additional resamplings that is done to balance the
% dataset, so that Nmur==NnoMur:
Nbalance = sum(~IposCases) - sum(IposCases);
% size of inflated data set (added samples to account for few positive samples):
NtrainTot = Ntrain + Nbalance;
clearvars Xtrain Ytrain
Xtrain = cell(NtrainTot,NsegDesired);
Ytrain = zeros(NtrainTot,NsegDesired);
% Mtrain = cell(NtrainTot,1);

for i=1:NtrainTot
    i 
    if i<=Ntrain
        ii = i;
        id = HSdataTrain.UNIKT_LOPENR(ii);
    else
        % randomly pick a sample of category 1 (done to balance data set)
        ii = randsample(JposCases,1);
        id = HSdataTrain.UNIKT_LOPENR(ii);
    end

    % get time series:
    x0 = wav2TS(id,aa);
    x = schmidt_spike_removal(x0, fs);
    x  = downsample(x,Nds);
    
    % get scaleogram and extract overlapping segement indeces:
    M = getScaleogram(x);
%     Mtrain{i} = M;
%     M = Mtrain{i};
    M = M(vertSlice,:);
    
    Nc_available = numel(segLines{ii,aa})-1;
    Ncs_actual = min(Ncs,Nc_available);
    [L,~] = getSegments(segLines{ii,aa},Nos,Ncs_actual,NsegDesired);
    
    % cycle through the segments:
    for k=1:NsegDesired
        % normalize, downsize, and save segment k:
        Mk = M(:,L(k,1):L(k,2));
        % normalize
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
        MkNorm = imresize(MkNorm, Ncomp);
        Xtrain{i,k} = MkNorm;
        Ytrain(i,k) = IposCases(ii);
    end 
end

Ytrain = categorical(Ytrain);
% reshape to a "vector" form:
Ytrain = reshape(Ytrain,[NsegDesired*NtrainTot,1]);
Xtrain = reshape(Xtrain,[NsegDesired*NtrainTot,1]);
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
load segLinesVal
%%

Nval = height(HSdataVal);
NsegVal = NsegDesired;
Xval = cell(Nval,NsegVal);
Yval = zeros(Nval,NsegVal);
Mval = cell(Nval,1);
% index for cases where there is murmur:
IposCases = HSdataVal.Murmur_1_grade_ref_ny_T72>=mg;
JposCases = find(IposCases);

for i=1:Nval
    i 
    
    id = HSdataVal.UNIKT_LOPENR(i);
    % get time series:
    x0 = wav2TS(id,aa);
    x = schmidt_spike_removal(x0, fs);
    x  = downsample(x,Nds);
    
    % get scaleogram and extract overlapping segement indeces:
    M = getScaleogram(x);
%     M = Mval{i};
    M = M(vertSlice,:);
    
    Nc_available = numel(segLinesVal{i,aa})-1;
    Ncs_actual = min(Ncs,Nc_available);
    [L,~] = getSegments(segLinesVal{i,aa},Nos,Ncs_actual,NsegDesired);
    
    % cycle through the segments:
    for k=1:NsegDesired
        % normalize, downsize, and save segment k:
        Mk = M(:,L(k,1):L(k,2));
        % normalize
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
        MkNorm = imresize(MkNorm, Ncomp);
        Xval{i,k} = MkNorm;
        Yval(i,k) = IposCases(i);
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
YtotPred = (sum(mode(reshape(Ypred,[Nval,NsegVal]),2)==mode(reshape(Yval,[Nval,NsegVal]),2)))/Nval
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
% 1 LSTM layer, 5 cps:[0.5385,0.9447]
% 1 LSTM layer, 5 cps:[0.6410,0.9112]
% 2 LSTM layer, 4 cps:[0.6154,0.9313]
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
clf
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

% using segmented data and downscaled scaleograms, with 3 segments per
% cycle, the algorithm achieves a validation score of 81% sensitivity and
% 60% specificity. This is now starting to look like some of the results
% that I saw in the physionet challenge. I used grade 2 as murmur
% threshold.

% I try to instead use 4 cycles per segment. I now get a sensitivity of 82%
% and a specificity of 72%; a very significant improvement. It appears that
% the LSTM algorithm BENEFITS from a recieving LONGER SEGMENTS as input.

% using 5 cycles per segment I achieved a sensitivity of 69% and a
% specificity of 85%. During training it seemed as if the algorithm was
% learning a lot faster, than with 4 cycles per segment, and so I thought
% the results would be better. It certainly achieved a higher training
% accuracy. Maybe I OVERTRAINED it? I will run the algorithm one more time,
% and try stopping it earlier, perhaps when it reaches 95% training
% accuracy.

% 59% sensitivity and 74% specificity was achieved when training was
% stopped earlier; clearly NOT an IMPROVEMENT. Stopping the algorithm
% earlier during training was not a good idea.
