% Cross validation and network training using only one position is used for
% training. Any type of target variable is accepted.

%#ok<*SEPEX>
%#ok<*NOPTS>
load nodes

Nsplits = 8;
% ¤¤ CHOOSE AUSCULTATION LOCATION ¤¤
aa = 2;
motherData = HSdata(union(Jtrain0,Jval0),:);
Iclean = motherData.(sprintf('MURMUR_%gNOISE_REF_T72',aa))==0;
% ¤¤ CHOOSE PRIMATRY OUTPUT VARIABLE THAT NETWORK PREDICTS ¤¤
% outputStr = sprintf('ARGRADE_T72');
outputStr = sprintf('Murmur_1_grade_ref_ny_T72');
% outputStr = "AVMEANPG_T72";
Javmean  = ~isnan(motherData.(outputStr));
Jclean       = find(and(Javmean,Iclean));
data0       = motherData(Jclean,:);
nodes0.loc      = nodes.loc(Jclean,:);
nodes0.state    = nodes.state(Jclean,:);
nodes0.seglines = nodes.seglines(Jclean,:);

% ¤¤ SHOULD TRAINING OUTPUT VARIABLE BE GRADED OR BINARY ¤¤
gradedOutput = true;
binaryOutput = ~gradedOutput;
if gradedOutput
    Y0 = data0.(outputStr);
elseif binaryOutput
    % ¤¤ IF NEEDED: CHOOSE POSITIVE CLASS THRESHOLD ¤¤
    classThr = 2;
    Y0 = data0.(outputStr)>=classThr;
end

N0   = height(data0);
Nval = floor(N0/Nsplits);
Jall = 1:N0;

% ¤¤ CHOOSE WHETHER OR NOT TO TRAIN NETWORK ¤¤
trainNet = true; % if false then gets validation info only
if trainNet
    networks = cell(1,Nsplits);
else not(trainNet)
    % ¤¤ IF NOT TRAINING: CHOOSE WHETHER OR NOT TO LOAD NETWORK ¤¤
    loadNetwork = true;
    if loadNetwork
    % ¤¤ IF LOAD NET: CHOOSE WHICH NETWORK TO LOAD ¤¤
%     load networksCVnoNoiseDrOutRegaa1234.mat
%     load networksCVnoNoiseDrOutaa4MRgeq1.mat
%     load networksCVnoNoiseDrOutARgeq1aa4.mat
        load networksCVnoNoiseDrOutaa1Reg.mat
    end
end

showPerfPlot = true;

AUCmat    = zeros(Nsplits,2);
NposCases = zeros(Nsplits,1);
CVresults.train.activations = cell(8,1);
CVresults.val.activations   = cell(8,1);
CVresults.train.J = cell(8,1);
CVresults.val.J   = cell(8,1);
Jpred  = zeros(Nval*Nsplits,1);

for i=1:Nsplits
    % ¤¤ SHOULD TRAINING BE INITIALIZED WITH PRETRAINED NETWORKS? ¤¤
    initTrainingWithNets = true;
    if initTrainingWithNets && trainNet
        initNets = load('networksCVnoNoiseDrOutMurRegaaAll.mat');
        initNet = initNets.networks{i};
    end
    
    i 
    % define the labels
    Jval   = (i-1)*Nval+1:i*Nval;
    Jtrain = setdiff(Jall,Jval);

    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
    
    % ¤¤ IF TRAINING: CHOOSE IF TO BALANCE DATA ¤¤
    balanceTrain = true;
    balanceVal   = true;
    Ncomp        = [30,200];

    if trainNet
        if ~gradedOutput
            posThr = 1;
        elseif gradedOutput
            % ¤¤ IF GRADED OUTPUT: SET THRESHOLD FOR POSITIVE CASES ¤¤
            posThr = 1; % all cases >= posThr are considered positive cases
        end
        
        [XtrainBal,YtrainBal,Xtrain,Ytrain] = genTrainOrValSet_new(data0,Y0,Jtrain,N,aa,...
                                    balanceTrain,nodes0,Ncomp,[],posThr);
    end
    
   [XvalBal,YvalBal,Xval,Yval] = genTrainOrValSet_new(data0,Y0,Jval,N,aa,...
                                    balanceVal,nodes0,Ncomp,[],posThr);
    
    % ¤¤ CHOOSE TRAINING TIME ¤¤
    trainTime = 60*60;
    miniBatchSize = 2^5;
    %%% training settings %%%
    if trainNet
        numFeatures = height(Xtrain{1}); % same as input size...
        inputSize = numFeatures; % number of timeseries to take as input
        numHiddenUnits = 50;
        if islogical(Y0)
            numClasses = 2;
            layers = [ ...
                sequenceInputLayer(numFeatures)
                lstmLayer(50)
                lstmLayer(50,'OutputMode','last')
                dropoutLayer(.5)
                fullyConnectedLayer(30)
                fullyConnectedLayer(numClasses)
                softmaxLayer
                classificationLayer];
        elseif isnumeric(Y0)
            layers = [ ...
                sequenceInputLayer(inputSize)
                lstmLayer(numHiddenUnits)
                lstmLayer(numHiddenUnits,'OutputMode','last')
                dropoutLayer(.5)
                fullyConnectedLayer(30)
                reluLayer
                fullyConnectedLayer(1)
                regressionLayer];
        end
        % ¤¤ CHOOSE MAX NUMBER OF EPOCHS:
        maxEpochs = 30;
        % ¤¤ SET INITIAL LEARN RATE AND LEARN-RATE-DROP-PERIOD ¤¤
        initLearnRate = 0.002;
        LearnRateDropPeriod = 5;
        validationFrequency = 498;%floor(numel(Ytrain)/miniBatchSize);
        % ¤¤ CHOOSE IF USE VAL-SET ACCURACY DURING TRAINING ¤¤
        checkValAccuracy = true;
        if checkValAccuracy
            validationData = {XvalBal,YvalBal};
        else
            validationData = [];
        end
        options = trainingOptions('adam', ...
            'ExecutionEnvironment','cpu', ...
            'MaxEpochs',maxEpochs, ...
            'MiniBatchSize',miniBatchSize, ...
            'GradientThreshold',1, ...
            'Plots','training-progress',...
            'LearnRateSchedule','piecewise', ...
            'LearnRateDropFactor',0.5, ...
            'ValidationData',validationData, ...
            'ValidationFrequency',validationFrequency, ...
            'LearnRateDropPeriod',LearnRateDropPeriod, ...
            'initialLearnRate',initLearnRate, ...
            'OutputFcn',@(info)myCostumOutputFcn(info,trainTime));
            % train network
%             load 'net_checkpoint__223__2021_12_18__19_52_29.mat' 'net';
            if initTrainingWithNets
                net = trainNetwork(XtrainBal,YtrainBal,initNet.Layers,options);
            else
                net = trainNetwork(XtrainBal,YtrainBal,layers,options);
            end
            % save network
            networks{i} = net;
    end
    
    Nloc = 1;
    if showPerfPlot
        figure
    end
    
    Ypred = predict(networks{i},Xtrain,'MiniBatchSize',miniBatchSize);
    CVresults.train.activations{i} = Ypred;
    
    Ypred = predict(networks{i},Xval,'MiniBatchSize',miniBatchSize);
    % ¤¤ CHOOSE VARIABLE TO PREDICT ¤¤
%     YpredTarget = modelData.Murmur_1_grade_ref_ny_T72(Jval)>=2;
    YpredTarget = data0.ASGRADE_T72(Jval)>=1;
    
    NposCases(i) = sum(YpredTarget);
%     Ypred2 = data0.Murmur_1_grade_ref_ny_T72(Jval)>=2;
    if NposCases(i)>0
        [AUC,X,Y,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                    size(Xval,1)/N.segPerPCG,N.segPerPCG,showPerfPlot,Nloc);
                    AUCmat(i,1) = AUC.seg;
                    AUCmat(i,2) = AUC.whole;
    else
        AUCmat(i,1) = nan;
        AUCmat(i,2) = nan;
        p.murWhole = median(reshapeSegmentScores(Ypred,numel(Jval),N.segPerPCG,Nloc),2);
    end
    CVresults.val.activations{i} = p.murWhole;
    CVresults.train.J{i} = Jclean(Jtrain);
    CVresults.val.J{i}   = Jclean(Jval);
    
    if showPerfPlot
        title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
        pause(.1)
    end

end

