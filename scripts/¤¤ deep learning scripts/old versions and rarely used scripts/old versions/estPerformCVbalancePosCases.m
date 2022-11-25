load murPredActivations_Reg_allAA.mat
% load murPredActivations_Reg_allAA.mat
%%
modelData = HSdata(union(Jtrain0,Jval0),:);
murValPredAll = paddedActivMatrix;
%% predictions using maximum murmur:
close all
Nsplits = 8;
% linear index for positions with atleast one clean signal to make
% prediction from:
% JatleastOne  = unique(cell2mat(S.joint.J));
% set the prediction set to be the set of auscultation areas where atleast
% one position has clean recording:
% predVar = 'ASGRADE_T72';
% predVar = 'Murmur_1_grade_ref_ny_T72';
predVar = 'maxMeanMurGrade';
classThr = 2;
predictionSet = modelData(JatleastOne,:);
maxMurActiv   = max(murValPredAll(JatleastOne,:),[],2);
JvalRedist = distributePosCasesEvenly(predictionSet.(predVar)>=classThr,...
                                                               Nsplits);
figures = false;
netName = 'netMurRegAllPos';
AUCmat.(netName).(predVar) = zeros(Nsplits,1);
for i=1:Nsplits
    Ytarget = predictionSet.(predVar)(JvalRedist.val{i})>=classThr;
    Yact    = maxMurActiv(JvalRedist.val{i});
    if figures; figure; end
    AUC = performanceSummaryNeurNet([],Ytarget,Yact,...
          numel(Yact),[],figures);
    if figures; title(sprintf('AUCwhole=%.3g, location=all',AUC.whole,aa)); end
    AUCmat.(netName).(predVar)(i) = AUC.whole;
end

computeCImeanEst(AUCmat.(netName).(predVar),'3','tDist')

%% get ROC curve for ALL max-murmur predictions:
% Prediction for entire data set
predVar = 'maxMeanMurGrade';
classThr = 2;
Ytarget = predictionSet.(predVar)>=classThr;
Yact    = maxMurActiv;
figure
[AUC] = performanceSummaryNeurNet([],Ytarget,Yact,...
    numel(Yact),[],true);
title(sprintf('AUCwhole=%.3g, location=%g',AUC.whole,aa));
%% get CV-AUC for each position
close all
Nsplits = 8;
for aa=1:4
% aa = 1;
% linear index for positions for which predictions have been made:
Jpred  = unique(cell2mat(S.indiv.J(:,aa)));
% set the prediction set to be the set of auscultation areas where atleast
% one position has clean recording:
predVar = 'ASGRADE_T72';
% predVar = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
classThr = 1;
predictionSet = modelData(Jpred,:);
maxMurActiv   = murValPredAll(Jpred,aa);
[JvalRedist.val,Npps] = distributePosCasesEvenly(predictionSet.(predVar)>=classThr,...
                                                               Nsplits);
figures = false;
netName = sprintf('netMurG2');
pos = sprintf('pos%g',aa);
AUCmat.(netName).(pos).(predVar).(num2name(classThr)) = zeros(Nsplits,1);

for i=1:Nsplits
    Ytarget = predictionSet.(predVar)(JvalRedist.val{i})>=classThr;
    Yact    = maxMurActiv(JvalRedist.val{i});
    if figures; figure; end
    AUC = performanceSummaryNeurNet([],Ytarget,Yact,...
        numel(Yact),[],figures);
    if figures; title(sprintf('AUCwhole=%.3g, location=%g',AUC.whole,aa));end
    AUCmat.(netName).(pos).(predVar).(num2name(classThr))(i) = AUC.whole;
end
end

% computeCImeanEst(AUCmat.(netName).(pos).(predVar).(num2name(classThr)),'3','tDist');
%% get ROC curve for all predictions combined, individual positions:
close all
aa = 1;
Jpred    = cell2mat(S.indiv.J(:,aa));
% predVar  = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
predVar  = 'ASGRADE_T72';
classThr = 1;
Ytarget = modelData.(predVar)(Jpred)>=classThr;
Yact    = cell2mat(S.indiv.activation(:,aa));
figure
[AUC,X,Y,T,p] = performanceSummaryNeurNet([],Ytarget,Yact,...
    numel(Yact),[],true);
title(sprintf('AUCwhole=%.3g, Ausc. Pos. = %g',AUC.whole,aa))
round(100*computeCImeanEst(AUC.whole,'2','tDist'),2)
%% test (p-values) G2 and regression nets against eachother
predVar  = 'ASGRADE_T72';
mean(AUCmat.netMurRegAllPos.(predVar) - AUCmat.netMurG2AllPos.(predVar))
pValue(AUCmat.netMurRegAllPos.(predVar) - AUCmat.netMurG2AllPos.(predVar))

predVar  = 'ARGRADE_T72';
xx = AUCmat.netMurRegPos4.(predVar) - AUCmat.netMurG2Pos4.(predVar);
mean(xx)
pValue(xx)

pVals    = zeros(3,4);
meanValDiffs = zeros(3,4);
%% Compute p-values for different class-thresholds and positions:
for aa=1:4
for classThr=1:3
% aa = 1;
% classThr = 1;
netName = sprintf('netMurG2');
pos = sprintf('pos%g',aa);
predVar = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
% predVar = 'ASGRADE_T72';
xx = AUCmat.netMurReg.(pos).(predVar).(num2name(classThr)) - ...
     AUCmat.netMurG2.(pos).(predVar).(num2name(classThr)) ;
meanValDiffs(classThr,aa) = mean(xx);
pVals(classThr,aa) = pValue(xx);
end
end


