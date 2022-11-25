% In this script, I estimate the optimal decision threshold based on each
% training set that was used during cross validation.
%% LOAD NETWORKS:
load networksCVnoNoiseDrOutRegaa1234.mat
load nodes
%% get training set performance and 
data0 = HSdata(union(Jtrain0,Jval0),:);
Jnonoise = data0.MURMUR_1NOISE_REF_T72==0;
data0 = data0(Jnonoise,:);
nodes0.loc      = nodes.loc(Jnonoise,:);
nodes0.state    = nodes.state(Jnonoise,:);
nodes0.seglines = nodes.seglines(Jnonoise,:);

% ¤¤ DEFINE LABELS FOR TRAINING DATA ¤¤
aa = 1;
mg = 2;
murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
Y0 = data0.(murSt)>=mg;

N0 = height(data0);
Nsplits = 8;
Nval = floor(N0/Nsplits);
Jall = 1:N0;

% ¤¤ CHOOSE IF GENERATE TRAINING DATA ¤¤
getTrainData = false;
if getTrainData
    activations = cell(Nsplits,1);
end
decisionVec = zeros(Nsplits,1);
AUCvecTrain = zeros(Nsplits,1);
AUCvecVal   = zeros(Nsplits,1);
sensVal = zeros(Nsplits,1);
specVal = zeros(Nsplits,1);
sensTrain = zeros(Nsplits,1);
specTrain = zeros(Nsplits,1);
uVec    = zeros(Nsplits,1);

% ¤¤ CHOOSE WHETHER OR NOT TO SHOW PERFORMANCE PLOTS ¤¤
figures = false;

for i=1:Nsplits
    if figures
        i %#ok<*NOPTS>
    end
    Jval   = (i-1)*Nval+1:i*Nval;
    Jtrain = setdiff(Jall,Jval);

    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;

    balanceTrain = false;
    balanceVal   = false;
    Ncomp = [30,200];
    
    if getTrainData
        [Xtrain,Ytrain] = genTrainOrValSet(data0,Y0,Jtrain,N,aa,...
                                 balanceTrain,nodes0,Ncomp,[]);
    end
    [Xval,Yval]     =  genTrainOrValSet(data0,Y0,Jval,N,aa,...
                                balanceVal,nodes0,Ncomp);
                            
    % get performance curve for training set and activations:
    Nloc = 1;
    Ntrain = size(Ytrain,1)/N.segPerPCG;
    if getTrainData
        miniBatchSize = 2^5;
        Ypred = predict(networks{i},Xtrain,'MiniBatchSize',miniBatchSize);
        activations{i} = Ypred;
    else
        Ypred = activations{i};
    end
    % ¤¤ CHOOSE VARIABLE TO PREDICT ¤¤
    YpredTarget = Y0(Jtrain);
    if figures; figure; end
    [AUC,X,Y,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                        Ntrain,N.segPerPCG,figures,Nloc);
                    
    minSensitivity = 0.93;
    u_opt = getOptimalThr(X,Y,T,minSensitivity);
    uVec(i) = -log(u_opt^-1-1);
    decisionVec(i) = u_opt;
    AUCvecTrain(i) = AUC.whole;
    Ypred = p.murWhole>=u_opt;
    sensTrain(i) = condProb(Ypred,YpredTarget);
    specTrain(i) = condProb(~Ypred,~YpredTarget);

    % compute corresponding sensitivity and specificity on validation set:
    
    Ypred = predict(networks{i},Xval,'MiniBatchSize',miniBatchSize);
    % ¤¤ CHOOSE VARIABLE TO PREDICT ¤¤
    YpredTarget = Y0(Jval);
    if figures; figure; end
    [AUC,X,Y,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
                    Nval,N.segPerPCG,figures,Nloc);
    AUCvecVal(i) = AUC.whole;
    % get predictions using above estimate of optimal threshold:
    Ypred = p.murWhole>=u_opt;
    % estimate sensitivity and specificity:
    sensVal(i) = condProb(Ypred,YpredTarget);
    specVal(i) = condProb(~Ypred,~YpredTarget);
    pause(.1)
end
