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
% number of cycles per segment:
Ncs = 4;
% cycle overlap for each consecutive pair of segments:
Nos = 2;
% number of segments to extract from each sample:
NsegPerPCG = 6;

Ncomp = [14,200];

XtrainAll = cell(4,1);
YtrainAll = cell(4,1);

for aa=1:4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% index for cases where there is murmur:
murStr = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = HSdataTrain.(murStr)>=mg;
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
        Ytrain(i,k) = IposCases(ii);
    end 
end

Ytrain = categorical(Ytrain);
% reshape to "vector" form:
Ytrain = reshape(Ytrain,[NsegPerPCG*NtrainTot,1]);
Xtrain = reshape(Xtrain,[NsegPerPCG*NtrainTot,1]);

% store the segments for auscultation location aa:
XtrainAll{aa} = Xtrain;
YtrainAll{aa} = Ytrain;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end
%% Combine into training segments into one set:
Xtrain = [XtrainAll{1};XtrainAll{2};XtrainAll{3};XtrainAll{4}];
Ytrain = [YtrainAll{1};YtrainAll{2};YtrainAll{3};YtrainAll{4}];
%% TEST SET
load segLinesVal
%%
Nval = height(HSdataVal);
XvalAll = cell(4,1);
YvalAll = cell(4,1);

for aa=1:4
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% index for cases where there is murmur:
murStr = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = HSdataVal.(murStr)>=mg;
JposCases = find(IposCases);
Xval = cell(Nval,NsegPerPCG);
Yval = zeros(Nval,NsegPerPCG);

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
        Mk = Mk(:,:);
        % normalize MFCC:
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
        MkNorm = imresize(MkNorm, Ncomp);
        Xval{i,k} = MkNorm;
        Yval(i,k) = IposCases(i);
    end 
end

% convert to type categorical:
Yval = categorical(Yval);
% reshape to a "vector" form:
Yval = reshape(Yval,[NsegPerPCG*Nval,1]);
Xval = reshape(Xval,[NsegPerPCG*Nval,1]);

% store the segments for auscultation location aa:
XvalAll{aa} = Xval;
YvalAll{aa} = Yval;
end
%% Combine validation segments from each location into one set:
Xval = [XvalAll{1};XvalAll{2};XvalAll{3};XvalAll{4}];
Yval = [YvalAll{1};YvalAll{2};YvalAll{3};YvalAll{4}];
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
    'initialLearnRate',0.002,...
    'OutputFcn',@(info)myCostumOutputFcn(info,7*3600));
%% Train network
net = trainNetwork(Xtrain,Ytrain,layers,options);
% save('trained networks/netAA1234murG2_v2.mat','net')
%% predict
Ypred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
figure
Nloc = 4;
[AUC,X,Y,T] = performanceSummaryNeurNet(Xval,Yval,net,Nval,NsegPerPCG,true,Nloc);

YpredTot = mode(reshapeSegmentScores(Ypred,Nval,NsegPerPCG,4),2);
YvalTot  = mode(reshapeSegmentScores(Yval, Nval,NsegPerPCG,4),2);
%% make predictions fpr a particular location:
k=1
sum(YpredTot((k-1)*Nval+1:k*Nval)==YvalTot((k-1)*Nval+1:k*Nval))/(1*Nval)
% sum(YpredWhole==YvalWhole)/(4*Nval)
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

%%%%% here begins comments from THIS script %%%%%%%%

% it turns out that despite being very slow to learn, and not appearing to
% achieve a good result when looking at training, using MFCC's has
% delivered the best performance yet, with a sensitivity and specificity of
% 82% and 77% respectively. It is possible that longer training would have
% improved the result.

% I now added a LSTM layer, for a total of 2 layers. This pushed max
% performance to 84.6% and 79.4% sensitivity and specificity respectively,
% the best results so far.

% I reduced the number of hidden units in the LSTM layers from 100 to 50.
% I got a significantly better looking ROC curve, with a peak result of
% 87.2% sensitivity and 77.4% specificity.

% I reduce the number of neurons in the last layer from 30 to 15. This
% resulted in an optimal performance of 79.5% sensitivity and 79.2%
% specificity, which is a significant REDUCTION in performance. Next I will
% try to increase the number of neurons instead.

% increasing the number of neurons in the fully connected layer did NOT
% IMPROVE classification performance. 30 neurons in the fully connected
% layer seems like a deceent choice.

% when using only 12 MFCC's I get what might be the best performance so
% far: 84.6% sensitivity and 80.1% specificity.



