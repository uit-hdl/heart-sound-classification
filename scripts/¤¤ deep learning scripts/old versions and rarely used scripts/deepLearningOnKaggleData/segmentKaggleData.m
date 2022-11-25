%% 
% load kaggle data:
ImportDataAndDefineVars
% load segmentation parameters:
loadHMMpar
%% GET LENGTH OF DATA
H = zeros(1,height(setb));
for i=1:height(setb)
    H(i) = numel(Xb{i});
end
% get index of recordings that are atleast 4 seconds long:
IsuffLong = H>3*4000;
%% GET EXTRACT DATA
clear XX
Idata = and(or(Ib.nor,Ib.mur),IsuffLong');
XX.all = Xb(Idata,1);
YY.all = setb.label(Idata)=="murmur";
%% SEGMENT DATA:
[~,Fs] = audioread(setb.fname(1));
Nds0 = 2;
Ntot = height(XX.all);

segLines = cell(Ntot,1);
for k=1:Ntot
    k %#ok<*NOPTS>
    x = XX.all{k};
    x = downsample(x,Nds0);
    fs = Fs/Nds0;
    assignedStates = runSpringerSegmentationAlgorithm(x, fs, HMMpar, [], false);
    [~,segLines{k}] = states2cycles(assignedStates);
 
%     clf
%     getScaleogram(x,1,true,false);
%     hold on
%     states2plot = [2,4];
%     cols        = ['r','g'];
%     title(sprintf('label:%s, sublabel:%s',setb.label(k),setb.sublabel(k)))
%     plotAssignedStates(assignedStates,states2plot,cols,1,50)
end
%% GET NUMBER OF CYCLES PER RECORDING

Ncycles = zeros(1,Ntot);
for i=1:Ntot
    Ncycles(i) = numel(segLines{i})-1;
end

XX.all   = XX.all(Ncycles>2);
YY.all   = YY.all(Ncycles>2);
segLines = segLines(Ncycles>2);

%% DIVIDE INTO TRAINING AND VALIDATION DATA
Ntot   = height(XX.all);
Ntrain = floor(Ntot*.8);
Nval  = Ntot - Ntrain;
Jtrain = randi(Ntot,1,Ntrain);
Jval = setdiff(1:Ntot,Jtrain);
segLinesTrain = segLines(Jtrain);
segLinesVal   = segLines(Jval);

XX.train = XX.all(Jtrain);
XX.val   = XX.all(Jval);
YY.train = YY.all(Jtrain);
YY.val   = YY.all(Jval);
%% CHECK THAT EXTRACTION OF TRAINING AND VALIDATION SETS WAS SUCCESSFUL:
k=1
clf
seg = segLines{k}(1):segLines{k}(2);
getScaleogram(XX.train{1}(seg),1,true,10:60,[1,inf]);
title(sprintf('%g',YY.val(k)))
%% generate training data and labels
fs = 4000;
% how much to downsample the signal:
Nds = 1;
% number of cycles per segment:
Ncs = 3;
% cycle overlap for each consecutive pair of segments:
Nos = 1;
% number of segments to extract from each sample:
NsegDesired = 3;
% to which size do we downsize the scaleogram segments:
Ncomp = [30,200];
% ySlice represents what vertical part of the scaleogram to use:
vertSlice = 10:53;

% index for cases where there is murmur:
IposCases = YY.all==1;
JposCases = find(IposCases);

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
    else
        % randomly pick a sample of category 1 (done to balance data set)
        ii = randsample(JposCases,1);
    end

    % get time series:
    x0 = XX.train{ii};
    x  = downsample(x0,Nds);
    
    % get scaleogram and extract overlapping segement indeces:
    M = getScaleogram(x,1,false,10:60,[1,inf]);
%     Mtrain{i} = M;
%     M = Mtrain{i};
    M = M(vertSlice,:);
    
    [L,~] = getSegments(segLinesTrain{ii},Nos,Ncs,NsegDesired);
    
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
Ytrain = reshape(Ytrain,[Ncs*NtrainTot,1]);
Xtrain = reshape(Xtrain,[Ncs*NtrainTot,1]);
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
NsegVal = NsegDesired;
Xval = cell(Nval,NsegVal);
Yval = zeros(Nval,NsegVal);
Mval = cell(Nval,1);
% index for cases where there is murmur:
IposCases = YY.train==1;
JposCases = find(IposCases);

for i=1:Nval
    i 

    % get time series:
    x0 = XX.val{i};
    x  = downsample(x0,Nds);
    
    % get scaleogram and extract overlapping segement indeces:
    M = getScaleogram(x);
%     M = Mval{i};
    M = M(vertSlice,:);
    [L,~] = getSegments(segLinesVal{i},Nos,Ncs,NsegVal);
    
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

% using segmented data and downscaled scaleograms, with 3 segments per
% cycle, the algorithm achieves a validation score of 81% sensitivity and
% 60% specificity. This is now starting to look like some of the results
% that I saw in the physionet challenge.


