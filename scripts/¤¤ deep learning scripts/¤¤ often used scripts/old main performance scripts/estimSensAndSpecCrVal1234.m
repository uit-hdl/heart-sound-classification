% Outdated script. Use CVscript... instead.

load nodes
motherData = HSdata(union(Jtrain0,Jval0),:);

Nsplits = 8;
Jall = union(Jtrain0,Jval0);
% ¤¤ CHOOSE NETWORK TO TEST ¤¤
load networksCVnoNoiseDrOutMurRegaaAll.mat

% ¤¤ CHOOSE WHETHER OR NOT TO GET TRAINING DATA ¤¤
getTrain = false;
if getTrain
    activTrainCell  = cell(Nsplits,4);
    YtrainCell      = cell(Nsplits,4);
else
    % ¤¤ IF ACTIVATIONS HAVE BEEN PREVIOUSLY TRAINED; LOAD THEM ¤¤
    loadActivations = false;
    if loadActivations
        load ..
    end
end
AUCmat.train  = zeros(Nsplits,4);
AUCmat.valAll = zeros(Nsplits,4);
AUCmat.valMax = zeros(Nsplits,1);
NposCases     = zeros(Nsplits,1);

murValPredAll = zeros(numel(Jall),4);

% define storage variables:
clear S
S.joint.J          = cell(8,1);
S.joint.activation = cell(8,1);
S.indiv.J          = cell(8,4);
S.indiv.activation = cell(8,4);
JtrainCell         = cell(8,4);
JvalCell           = cell(8,4);

% *** plotting options ***
plotTrain    = false;
plotValIndiv = false;
plotValJoint = false;

for i=1:Nsplits
i %#ok<*NOPTS>
Xtrain = cell(4,1);
Ytrain = cell(4,1);
Xval   = cell(4,1);
Yval   = cell(4,1);

for aa=1:4
    aa
    % set the mother-data set:
    data0  = motherData;
    Jclean = find(data0.(sprintf('MURMUR_%gNOISE_REF_T72',aa))==0);
    data0    = data0(Jclean,:);
    nodes0.loc      = nodes.loc(Jclean,:);
    nodes0.state    = nodes.state(Jclean,:);
    nodes0.seglines = nodes.seglines(Jclean,:);
    
    % define the labels
    murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    % ¤¤ IF TRAINING: CHOOSE THE OUTPUT DATA TO TRAIN ON ¤¤
    Y0   = data0.(murSt)>=2;
    N0   = height(data0);
    Nval = floor(N0/Nsplits);
    Jall = 1:N0;
    
    Jval   = (i-1)*Nval+1:i*Nval;
    Jtrain = setdiff(Jall,Jval);
    
    clear N
    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
    
    balanceTrain = false;
    balanceVal   = false;
    Ncomp = [30,200];
    
    if getTrain
       [Xtrain{aa},Ytrain{aa}] = genTrainOrValSet(data0,Y0,Jtrain,N,aa,...
                                        balanceTrain,nodes0,Ncomp);
    end
    
   [Xval{aa},Yval{aa}]     =  genTrainOrValSet(data0,Y0,Jval,N,aa,...
                                    balanceVal,nodes0,Ncomp);
    
    JtrainCell{i,aa} = Jclean(Jtrain);
    JvalCell{i,aa}   = Jclean(Jval);
end
% combine into 1 training set and validation set:
if getTrain
    Xtrain1234 = [Xtrain{1};Xtrain{2};Xtrain{3};Xtrain{4}];
    Ytrain1234 = [Ytrain{1};Ytrain{2};Ytrain{3};Ytrain{4}];
end
Xval1234   = [Xval{1};Xval{2};Xval{3};Xval{4}];
Yval1234   = [Yval{1};Yval{2};Yval{3};Yval{4}];

miniBatchSize = 2^5;

% *** training set predictions ***
if plotTrain
for aa=1:4
    if getTrain
        Ypred = predict(networks{i},Xtrain{aa},'MiniBatchSize',miniBatchSize);
    else
        Ypred = activTrainCell{i,aa};
    end
    % ¤¤ IF TRAINING: CHOOSE TRAINING SET PREDICITON TARGET ¤¤
    predVar = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    YpredTarget = motherData.(predVar)(JtrainCell{i,aa})>=2;
    subplot(2,2,aa)
    [AUC,X,Y,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        numel(Ypred(:,end))/N.segPerPCG,N.segPerPCG,true); %#ok<*ASGLU>
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
    AUCmat.train(i,aa) = AUC.whole;
    
    YtrainCell{i,aa} = YpredTarget;
    if getTrain
        % save training activations:
        activTrainCell{i,aa} = p.murWhole;
    end
end
pause(.5)
end

% *** validation predictions for each location ***
if plotValIndiv
    figure
end
for aa=1:4
    Ypred = predict(networks{i},Xval{aa},'MiniBatchSize',miniBatchSize);
    % ¤¤ CHOOSE PREDICITON TARGET ¤¤ 
    predVar = sprintf('Murmur_%G_grade_ref_ny_T72',aa);
    classThr = 2;
    YpredTarget = motherData.(predVar)(JvalCell{i,aa})>=classThr;
    
    subplot(2,2,aa)
    [AUC,X,Y,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        numel(Ypred(:,end))/N.segPerPCG,N.segPerPCG,plotValIndiv);
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
    
    AUCmat.valAll(i,aa) = AUC.whole;
    murValPredAll(JvalCell{i,aa},aa) = p.murWhole;
    
    % save prediction info for validation set i, location aa:
    S.indiv.J{i,aa}          = JvalCell{i,aa};
    S.indiv.activation{i,aa} = p.murWhole;
end
pause(.2)

% *** max murmur predictions ***
% YpredTarget = dataModel.maxMeanMurGrade(Jval)>=2;
JatleastOne = unionIterated(JvalCell(i,:));
predVar = 'ARGRADE_T72';
classThr = 3;
YpredTarget = motherData.(predVar)(JatleastOne)>=classThr;
Ypred       = max(murValPredAll(JatleastOne,:),[],2);
NposCases(i) = sum(YpredTarget);
if NposCases(i)>0
    if plotValJoint
        figure
    end
    [AUC,X,Y,T,p] = performanceSummaryNeurNet([],YpredTarget,Ypred,...
        numel(Ypred(:,end)),N.segPerPCG,plotValJoint);
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=all',AUC.seg,AUC.whole,aa))
    AUCmat.valMax(i) = AUC.whole;
else
    AUCmat.valMax(i) = nan;
end
pause(.2)

% save adress and activations for max predictions:
S.joint.J{i}          = JatleastOne;
S.joint.activation{i} = Ypred;

clearvars Xtrain1234 Ytrain1234 Xval1234 Yval1234

end

%
% CVresults.val.activations    = S.indiv.activation;
% CVresults.val.J              = JvalCell;
% CVresults.train.activations  = activTrainCell;
% CVresults.train.J            = JtrainCell;
% save CVresults_netMurG2AllPos.mat CVresults
