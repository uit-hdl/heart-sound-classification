% In this script, I estimate the optimal decision threshold based on each
% training set that was used during cross validation.
%% LOAD NETWORKS:
load networksCVnoNoiseDrOutReg.mat
load nodes
%% get training set performance and 
data0    = HSdataTrainAndVal;
Jnonoise = data0.MURMUR_1NOISE_REF_T72==0;
data0    = data0(Jnonoise,:);
nodes0.loc      = nodes.loc(Jnonoise,:);
nodes0.state    = nodes.state(Jnonoise,:);
nodes0.seglines = nodes.seglines(Jnonoise,:);

% define the labels
aa = 1;
mg = 2;
murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = data0.(murSt)>=mg;
N0 = height(data0);
Nsplits = 8;
Nval = floor(N0/Nsplits);
Jall = 1:N0;
activations = cell(Nsplits,1);
decisionVec = zeros(Nsplits,1);
AUCvecTrain = zeros(Nsplits,1);
AUCvecVal   = zeros(Nsplits,1);
sensVec = zeros(Nsplits,1);
specVec = zeros(Nsplits,1);
uVec    = zeros(Nsplits,1);

plotTrain = false;
plotVal   = false;

for i=1:Nsplits
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

    YgradeTrain = Ygrade(Jtrain);
    YgradeVal   = Ygrade(Jval);
    balanceTrain = true;
    balanceVal   = false;
    Ncomp = [30,200];

    % get training set and validation set
    [Xtrain,Ytrain] = genTrainOrValSet(dataTrain,N,aa,YgradeTrain,...
                            balanceTrain,segLinesTrain,Ncomp);
    [Xval,Yval] = genTrainOrValSet(dataVal,N,aa,YgradeVal,balanceVal,...
                            segLinesVal,Ncomp);
%     reshape(Yval,[Nval,N.segPerPCG]);
%     mode(reshapeSegmentScores(Yval>=mg,Nval,N.segPerPCG,Nloc),2)

    % get performance curve for training set and activations:
    Nloc = 1;
    Ntrain = size(Ytrain,1)/N.segPerPCG;
    miniBatchSize = 2^5;
    YpredTrain = predict(networks{i},Xtrain,'MiniBatchSize',miniBatchSize);
    figure
    [AUC,X,Y,T,p] = performanceSummaryNeurNet(Xtrain,Ytrain>=mg,YpredTrain,...
                            Ntrain,N.segPerPCG,plotTrain,Nloc);
    minSensitivity = 0.9;
    u_opt = getOptimalThr(X,Y,T,minSensitivity);
    % apply inverse logistic function to convert threshold to murmur grade.
%     u_opt = -log(u_opt^-1-1);
    uVec(i) = u_opt;
    decisionVec(i) = u_opt;
    activations{i} = p.murWhole;
    AUCvecTrain(i) = AUC.whole;

    % compute corresponding sensitivity and specificity on validation set:
    YpredReg = predict(networks{i},Xval,'MiniBatchSize',miniBatchSize);
    figure
    [AUC,X,Y,T,p,YvalWhole] = performanceSummaryNeurNet(Xval,Yval>=mg,YpredReg,...
                                        numel(Jval),N.segPerPCG,plotVal,Nloc);
     getOptimalThr(X,Y,T,minSensitivity);
     AUCvecVal(i) = AUC.whole;
     % get predictions using above estimate of optimal threshold:
     YpredVal = p.murWhole>=u_opt;
     % estimate sensitivity and specificity:
     sensVec(i) = condProb(YpredVal,YvalWhole)
     specVec(i) = condProb(~YpredVal,~YvalWhole)
     uVec
     pause(.1)
end
%%
uVec = -log(uVec.^-1 - 1);
save('performanceCVnoNoiseDropOutReg',...
    'AUCvecTrain','AUCvecVal','sensVec','specVec','uVec')
