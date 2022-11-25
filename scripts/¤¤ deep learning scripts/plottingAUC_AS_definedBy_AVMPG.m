%% estimate PREDICTIVE power of MAX-MURMUR-PREDICTION algorithm
% ¤¤ SET NUMBER OF CV-SPLITS ¤¤
net = "regOver";
if net=="reg"
    load CVresults_netMurRegAllPos_valStop
elseif net=="regOver"
    load CVresults_netMurRegAllPos_valStop_overTrain
elseif net=="G2"
    load CVresults_netMurG2AllPos_valStop
elseif net=="G2Over"
    load CVresults_netMurG2AllPos_valStop_overTrain
end

%% Preliminary:
NrowsTot = height(HSdata);
predMatrix   = zeros(NrowsTot,4);
targetMatrix = zeros(NrowsTot,4);
sympAndSickMat = zeros(NrowsTot,4);
noSympAndSickMat = zeros(NrowsTot,4);

%%
Nsplits = 8;
XX = zeros(3,21);
for k=10:30
% define storage variables:
sn    = zeros(Nsplits,1);
sp    = zeros(Nsplits,1);
ac    = zeros(Nsplits,1);
AUCmat = zeros(Nsplits,1);
predictions = cell(Nsplits,1);
activations = cell(Nsplits,1);
Ytrue = cell(Nsplits,1);
Jall = cell2mat(CVresults.valTot.J);

% allVHDpredMatrix = zeros(height(HSdata),4);
% plot settings:
plotTrain = false;
plotVal   = false;
plotValTot = true;
plotComparisonSubplot = false;
if ~plotComparisonSubplot
    close all
end
for i=1:Nsplits
    targetType = 'AS';
    HSdata.ASPGgrade = HSdata.AVMEANPG_T72>=k;
    if targetType(2)=='R'
        classThr = 3;
    else
        classThr = 1;
    end
    targetVar = sprintf('%sgrade',targetType);
    I_disease = disease2index(targetType);
    
    ActMatVal = getZeroPaddedActivMatrix(CVresults.val.activations(i,:),...
                                   CVresults.val.J(i,:),NrowsTot);
    ActMatTrain = getZeroPaddedActivMatrix(CVresults.train.activations(i,:),...
                                   CVresults.train.J(i,:),NrowsTot);

    plotVal = false;
    ActMat = ActMatTrain + ActMatVal;
    minSn = 0.5;
    minSp = 0.5;
    [activ,u0,Ytarget,Ypred,AUC] = get_sigVHDactivations(ActMat,HSdata,...
                                            CVresults.trainTot.I{i},...
                                            CVresults.valTot.I{i},...
                                            targetType,classThr,...
                                            plotVal,minSn,minSp);
    % store results:                               
    AUCmat(i)      = AUC;
    predictions{i} = Ypred.val;
    activations{i} = activ.val;
    Ytrue{i}       = Ytarget.val;
    
    ac(i) = mean((activ.val>=u0)==Ytarget.val);
    sn(i)  = condProb(Ypred.val,Ytarget.val);
    sp(i)  = condProb(~Ypred.val,~Ytarget.val);
    
    targetMatrix(CVresults.valTot.I{i},I_disease) = Ytarget.val;
    predMatrix(CVresults.valTot.I{i},I_disease)   = Ypred.val;
    
end

predictions = cell2mat(predictions);
YpredPadded = predMatrix(:,I_disease);
activations = cell2mat(activations);
YtargetTot = HSdata.avmeanpg(Jall)>=k;

if plotValTot
    hold on
    AUCtot = performanceSummaryNeurNet([],YtargetTot,activations,...
                              [],[],plotValTot);
end
ci = getBootStrappedAUCci(YtargetTot,activations,100);
% myFunc = @(x1,x2) log(getAUC(x1,x2));
% BOOTSTAT = bootstrp(200,myFunc,YtargetTot,activations);
% bootStd = std(BOOTSTAT);
% ci = log(AUCtot.whole) + [-1 1]*bootStd*1.96;
% ci = exp(ci);
XX(:,k-9) = [ci(1);AUCtot.whole;ci(2)]
% k
% myFunc = @(x1,x2) getAUC(x1,x2);
% myFunc = @(x1,x2) condProb(x2>=.3,x1);
% nIterations = 100;
% compute the CI of correlation between murmur grade and
% aortic-valve-pressure-gradient:
% CI = bootci(nIterations,...
%     {myFunc, YtargetTot, activations})
end
%%
close
plotWithCI(100*(XX),10:30,["grey","k","grey"],[],["95% boostrap-ci","AUC"]);
xlabel 'AVMPG-threshold (mm Hg)'
ylabel 'AUC (%)'



