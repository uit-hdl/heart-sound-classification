% generate validation and training sets:
load nodes


NvalSets = 8;
networks = cell(1,8);

for i=1:NvalSets
    
Xtrain = cell(4,1);
Ytrain = cell(4,1);
Xval = cell(4,1);
Yval = cell(4,1);

for aa=1:4
    data0    = HSdataTrainAndVal;
    Jnonoise = data0.(sprintf('MURMUR_%gNOISE_REF_T72',aa))==0;
    data0    = data0(Jnonoise,:);
    nodes0.loc      = nodes.loc(Jnonoise,:);
    nodes0.state    = nodes.state(Jnonoise,:);
    nodes0.seglines = nodes.seglines(Jnonoise,:);

    % define the labels
    murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    Ygrade = data0.(murSt);
    N0 = height(data0);
    Nval = floor(N0/NvalSets);
    Jall = 1:N0;


    i %#ok<*NOPTS>
    Jval   = (i-1)*Nval+1:i*Nval;
    Jtrain = setdiff(Jall,Jval);
    dataTrain = data0(Jtrain,:);
    dataVal   = data0(Jval  ,:);
    segLinesTrain = nodes0.seglines(Jtrain,:);
    segLinesVal   = nodes0.seglines(Jval  ,:);

    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;

    YgradeTrain = Ygrade(Jtrain);
    YgradeVal   = Ygrade(Jval);
    balanceTrain = true;
    balanceVal   = false;
    Ncomp = [30,200];

   [Xtrain{aa},Ytrain{aa}] = genTrainOrValSet(dataTrain,N,aa,YgradeTrain,...
                                balanceTrain,segLinesTrain,Ncomp);
   [Xval{aa},Yval{aa}] = genTrainOrValSet(dataVal,N,aa,YgradeVal,balanceVal,...
                                segLinesVal,Ncomp);
end
% combine into 1 training set and validation set:
Xtrain1234 = [Xtrain{1};Xtrain{2};Xtrain{3};Xtrain{4}];
Ytrain1234 = [Ytrain{1};Ytrain{2};Ytrain{3};Ytrain{4}];
Xval1234 = [Xval{1};Xval{2};Xval{3};Xval{4}];
Yval1234 = [Yval{1};Yval{2};Yval{3};Yval{4}];

%%% train %%%
numFeatures = height(Xtrain{1}{1}); % same as input size...
inputSize = numFeatures; % number of timeseries to take as input
numHiddenUnits = 50;
numClasses = 2;
layers = [ ...
sequenceInputLayer(inputSize)
lstmLayer(numHiddenUnits)
lstmLayer(numHiddenUnits,'OutputMode','last')
dropoutLayer(.5)
fullyConnectedLayer(30)
reluLayer
fullyConnectedLayer(1)
regressionLayer];
% specify training options
maxEpochs = 70;
miniBatchSize = 2^5;
options = trainingOptions('adam', ...
    'ExecutionEnvironment','cpu', ...
    'MaxEpochs',maxEpochs, ...
    'MiniBatchSize',miniBatchSize, ...
    'GradientThreshold',1, ...
    'Plots','training-progress',...
    'LearnRateSchedule','piecewise', ...
    'LearnRateDropFactor',0.5, ...
    'LearnRateDropPeriod',5, ...
    'initialLearnRate',0.002, ...   
    'OutputFcn',@(info)myCostumOutputFcn(info,270*60,inf));
% train network
net = trainNetwork(Xtrain1234,Ytrain1234,layers,options);
% save network
networks{i} = net;

Nloc = 1;
figure
for aa=1:4
    Ypred = predict(net,Xval{aa},'MiniBatchSize',miniBatchSize);
    subplot(2,2,aa)
    [AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval{aa},Yval{aa}>=2,Ypred,...
        size(Xval{aa},1)/N.segPerPCG,N.segPerPCG,true,Nloc);
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
    pause(.1)
end
clearvars Xtrain1234 Ytrain1234 Xval1234 Yval1234
end

%% ESTIMATE SENSITIVITY AND SPECIFICITY
% produce validation set
for i=1:8
    i
%     i = 1; %#ok<*NOPTS>
    Jval   = (i-1)*Nval+1:i*Nval;
    dataVal   = data0(Jval,:);
    segLinesVal = nodes0.seglines(Jval,:);
    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
    YgradeVal = Ygrade(Jval);
    balanceVal   = false;
    Ncomp = [30,200];
    [Xval,Yval] = genTrainOrValSet(dataVal,N,aa,YgradeVal,balanceVal,...
                                    segLinesVal,Ncomp);

    [AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval,Yval,networks{i},numel(Jval),N.segPerPCG,false,Nloc);

    NpcgPerLoc = Nval;
    NsegPerPcg = 6;
    Nloc = 1;

    opt_x_coord = .9; %get x-coordinate of a good point on performance curve
    [~,imin] = min(abs(Y.whole - opt_x_coord))
%     u = T.whole(imin)
    u = 0.65
    % Make predictions with new prediction threshold:
    YpredTot = categorical(double(p.murWhole > u*(1 - .1)));
    YvalTot  = mode(reshapeSegmentScores(Yval,Nval,NsegPerPcg,Nloc),2);
    % compute accuracy:
    sum(YpredTot==YvalTot)/Nval
    % sensitivity and specificity:
    ss(i,:) = [condProb(YpredTot=='1',YvalTot=='1'),...
    condProb(YpredTot=='0',YvalTot=='0')] %#ok<SAGROW>
end
%% whole segments:
dataVal = HSdataVal;
clearvars I N
mur   = categorical(1);
nomur = categorical(0);
I.pred.mur   = (YpredTot==mur)';
I.pred.nomur = (YpredTot==nomur)';
I.actual.mur = (YvalTot==mur)';
I.actual.nomur = (YvalTot==nomur)';

I.truePos  = and(I.pred.mur,I.actual.mur);   % correctly classified murmur
I.falsePos = and(I.pred.mur,I.actual.nomur); % incorrectly classified as murmur
I.trueNeg  = and(I.pred.nomur,I.actual.nomur);
I.falseNeg = and(I.pred.nomur,I.actual.mur); % murmurs that were missed
N.tot      = height(dataVal);
N.mur      = sum(I.actual.mur);
N.nomur      = sum(I.actual.nomur);
N.falsePos = sum(I.falsePos);
N.falseNeg = sum(I.falseNeg);
N.truePos  = sum(I.truePos); %#ok<*NOPTS>
N.trueNeg  = sum(I.trueNeg)
accuracy = (N.truePos + N.trueNeg)/N.tot

