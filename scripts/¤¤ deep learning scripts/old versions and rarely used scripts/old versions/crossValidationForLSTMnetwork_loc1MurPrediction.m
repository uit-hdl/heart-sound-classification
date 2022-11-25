% generate validation and training sets:
load nodes

data0 = HSdataTrainAndVal;
Jnonoise = data0.MURMUR_1NOISE_REF_T72==0;
data0 = data0(Jnonoise,:);
nodes0.loc      = nodes.loc(Jnonoise,:);
nodes0.state    = nodes.state(Jnonoise,:);
nodes0.seglines = nodes.seglines(Jnonoise,:);

% define the labels
aa = 1;
mg = 2;
murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = data0.(murSt)>=mg;
N0 = height(data0);
NvalSets = 8;
Nval = floor(N0/NvalSets);
Jall = 1:N0;

getPerformance = false;
% cell array that stores networks
if getPerformance
    load networks
else
    networks = cell(1,8);
end

for i=1:NvalSets
i %#ok<*NOPTS>
Jval   = (i-1)*Nval+1:i*Nval;
Jtrain = setdiff(Jall,Jval);
dataTrain = data0(Jtrain,:);
dataVal   = data0(Jval,:);
segLinesTrain = nodes0.seglines(Jtrain,:);
segLinesVal   = nodes0.seglines(Jval,:);

N.ds = 20;
N.cs = 4;
N.os = 2;
N.segPerPCG = 6;

IposCasesTrain = IposCases(Jtrain);
IposCasesVal   = IposCases(Jval);
balanceTrain = true;
balanceVal   = false;
Ncomp = [30,200];


if ~getPerformance
   [Xtrain,Ytrain] = genTrainOrValSet(dataTrain,N,aa,IposCasesTrain,...
                                balanceTrain,segLinesTrain,Ncomp);
   [Xval,Yval] = genTrainOrValSet(dataVal,N,aa,IposCasesVal,balanceVal,...
                                segLinesVal,Ncomp);
else
    [Xval,Yval] = genTrainOrValSet(dataVal,N,aa,IposCasesVal,balanceVal,...
                                segLinesVal,Ncomp);
end

if ~getPerformance
%%% train %%%
numFeatures = height(Xtrain{1}); % same as input size...
numHiddenUnits = 50;
numClasses = 2;
layers = [ ...
    sequenceInputLayer(numFeatures)
    lstmLayer(numHiddenUnits)
    lstmLayer(numHiddenUnits,'OutputMode','last')
    dropoutLayer(.5)
    fullyConnectedLayer(30)
    fullyConnectedLayer(numClasses)
    softmaxLayer
    classificationLayer];
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
        'LearnRateDropPeriod',5,...
        'initialLearnRate',0.002, ...   
        'OutputFcn',@(info)myCostumOutputFcn(info,70*60,30));
    
    % train network
    net = trainNetwork(Xtrain,Ytrain,layers,options);
    % save network
    networks{i} = net;
else
    net = networks{i};
end

Nloc = 1;
figure()
[AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval,Yval,net,numel(Jval),N.segPerPCG,true,Nloc);
title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
pause(.1)
end

%% ESTIMATE SENSITIVITY AND SPECIFICITY
% produce validation set
for i=1:8
    i
%     i = 1; %#ok<*NOPTS>
    Jval      = (i-1)*Nval+1:i*Nval;
    dataVal   = data0(Jval,:);
    segLinesVal = nodes0.seglines(Jval,:);
    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
    IposCasesVal = IposCases(Jval);
    balanceVal   = false;
    Ncomp = [30,200];
    [Xval,Yval] = genTrainOrValSet(dataVal,N,aa,IposCasesVal,balanceVal,...
                                    segLinesVal,Ncomp);

    [AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval,Yval,networks{i},...
                                numel(Jval),N.segPerPCG,false,Nloc);

    NpcgPerLoc = Nval;
    NsegPerPcg = 6;
    Nloc = 1;

    opt_x_coord = .9; %get x-coordinate of a good point on performance curve
    [~,imin] = min(abs(Y.whole - opt_x_coord))
%     u = T.whole(imin)
    u = 0.65
    % Make predictions with new prediction threshold:
    YpredTot = categorical(double(p.murWhole > u*(1 - 0)));
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
N.nomur    = sum(I.actual.nomur);
N.falsePos = sum(I.falsePos);
N.falseNeg = sum(I.falseNeg);
N.truePos  = sum(I.truePos); %#ok<*NOPTS>
N.trueNeg  = sum(I.trueNeg)
accuracy = (N.truePos + N.trueNeg)/N.tot

