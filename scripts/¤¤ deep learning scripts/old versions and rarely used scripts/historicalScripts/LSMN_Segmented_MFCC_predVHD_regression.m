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
NsegPerPCG = 12;
L = getSegments(segLines{6,1},Nos,Ncs,NsegPerPCG);
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
aa = 4;
% number of cycles per segment:
Ncs = 4;
% cycle overlap for each consecutive pair of segments:
Nos = 2;
% number of segments to extract from each sample:
NsegPerPCG = 6;
% to which size do we downsize the scaleogram segments:
Ncomp = [30,200];
% ySlice represents what vertical part of the scaleogram to use:
vertSlice = 10:53;  

% index for cases where there is murmur:
vhd = 'MRGRADE_T72';
grades = HSdataTrain.(vhd);
IposCases = grades>0;
JposCases = find(IposCases);

Ntrain = height(HSdataTrain);
% Nbalance is number of additional resamplings that is done to balance the
% dataset, so that Nmur==NnoMur:
Nbalance = sum(~IposCases) - sum(IposCases);
% size of inflated data set (added samples to account for few positive samples):
NtrainTot = Ntrain + Nbalance;
clearvars Xtrain Ytrain
Xtrain = cell(NtrainTot,NsegPerPCG);
Ytrain = zeros(NtrainTot,NsegPerPCG);
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
    
    Nc_available = numel(segLines{ii,aa})-1;
    Ncs_actual = min(Ncs,Nc_available);
    [L,~] = getSegments(segLines{ii,aa},Nos,Ncs_actual,NsegPerPCG);
    
    % cycle through the segments:
    for k=1:NsegPerPCG
        % get the segment:
        xk = x(L(k,1):L(k,2));
        % get compact representation of signal in time-frequency domain:
        Mk = getMFCC(xk,floor(fs/Nds));
        % extract the number of coefficients desired as features:
        Mk = Mk(:,:);
        % normalize MFCC:
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
        MkNorm = imresize(MkNorm, Ncomp);
        Xtrain{i,k} = MkNorm;
        Ytrain(i,k) = grades(ii);
    end 
end

% reshape to a "vector" form:
Ytrain = reshape(Ytrain,[NsegPerPCG*NtrainTot,1]);
Xtrain = reshape(Xtrain,[NsegPerPCG*NtrainTot,1]);
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
Xval = cell(Nval,NsegPerPCG);
Yval = zeros(Nval,NsegPerPCG);
% index for cases where there is murmur:
grades = HSdataVal.(vhd);
IposCases = grades>0;
JposCases = find(IposCases);

for i=1:Nval
    i 
    
    id = HSdataVal.UNIKT_LOPENR(i);
    % get time series:
    x0 = wav2TS(id,aa);
    x = schmidt_spike_removal(x0, fs);
    x  = downsample(x,Nds);
    
    Nc_available = numel(segLinesVal{i,aa})-1;
    Ncs_actual = min(Ncs,Nc_available);
    [L,~] = getSegments(segLinesVal{i,aa},Nos,Ncs_actual,NsegPerPCG);
    
    % cycle through the segments:
    for k=1:NsegPerPCG
        % get the segment:
        xk = x(L(k,1):L(k,2));
        % get compact representation of signal in time-frequency domain:
        Mk = getMFCC(xk,floor(fs/Nds));
        % extract the number of coefficients desired as features:
        Mk = Mk(1:14,:);
        % normalize MFCC:
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
        MkNorm = imresize(MkNorm, Ncomp);
        Xval{i,k} = MkNorm;
        Yval(i,k) = IposCases(i);
    end 
end

Yval = double(Yval);
% reshape to a "vector" form:
Yval = reshape(Yval,[NsegPerPCG*Nval,1]);
Xval = reshape(Xval,[NsegPerPCG*Nval,1]);
%% Define LSTM architecture
numFeatures = height(Xtrain{1}); % same as input size...
inputSize = numFeatures; % number of timeseries to take as input

numHiddenUnits = 50;
numClasses = 2;

layers = [ ...
    sequenceInputLayer(inputSize)
    lstmLayer(numHiddenUnits)
    lstmLayer(numHiddenUnits,'OutputMode','last')
    fullyConnectedLayer(30)
    reluLayer
    fullyConnectedLayer(1)
    regressionLayer];

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
    'initialLearnRate',0.002);
%% Train network
net = trainNetwork(Xtrain,Ytrain,layers,options);
%% predict
Ypred = classify(net,Xval,'MiniBatchSize',miniBatchSize);

Ypred = predict(net,Xval,'MiniBatchSize',miniBatchSize)

% take most common prediction of data-sample as the prediction
YpredTot = mode(reshape(Ypred,[Nval,NsegPerPCG]),2);
YvalTot  = mode(reshape(Yval,[Nval,NsegPerPCG]),2);

sum(YpredTot==YvalTot)/Nval
%%
clear I
mur   = categorical(1);
noMur = categorical(0);
I.pred.mur   = (Ypred==mur)';
I.val.mur    = (Yval==mur)';

acc = sum(Ypred == Yval)./numel(Yval)
accTot = sum(YpredTot == YvalTot)./numel(YvalTot)

NcorrPredAndMur = sum(and(Ypred == Ypred,Ypred==mur));

sens = sum(and(Ypred == Yval,Ypred==mur))/sum(Yval==mur)
spec = sum(and(Ypred==Yval,Ypred==noMur))/sum(Yval==noMur)
%%
sensTot = sum(and(YpredTot == YvalTot,YpredTot==mur))/sum(YvalTot==mur)
specTot = sum(and(YpredTot==YvalTot,YpredTot==noMur))/sum(YvalTot==noMur)

% 1 LSTM layer, 5 cps:[0.5385,0.9447]
% 1 LSTM layer, 5 cps:[0.6410,0.9112]
% 2 LSTM layer, 4 cps:[0.6154,0.9313]
% 1 LSTM layer, 4 cps, MFCC(all):[0.5897,0.8878] 
% 2 LSTM layer, 4 cps, MFCC(all):[0.7949,0.8593] 
% 2 LSTM layer, 4 cps, MFCC(all), 50 hidden units:[0.7949,0.8559] 
% 2 LSTM layer, 4 cps, MFCC(all), 50 hidden units, 15 neurons in fully connected :[0.6667,0.8492]
% 2 LSTM layer, 4 cps, MFCC(all), 50 hidden units, 15 neurons in fully connected :[0.6923,0.8074] 
%% COMPUTE RUC CURVE
act = activations(net,Xval,'softmax');

I.pred.mur = act(2,:)>.01
xx = 0.00001:.0002:0.99;
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