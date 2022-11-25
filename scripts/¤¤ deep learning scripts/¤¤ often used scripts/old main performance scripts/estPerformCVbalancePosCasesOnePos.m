% outdated script. Computes performance metrics using cross validation
% results.

%% load the validation activations
% load murPredActivations_Reg_allAA.mat
% load murPredActivations_G2_allAA.mat
% load CVresults_netMurRegAllPos.mat
load CVresults_AVMPG_firstTry.mat
%% preliminary
modelData = HSdataTrain;
predictionSet  = modelData(Jclean,:);
activationsVal = cell2mat(CVresults.val.activations(:,1));
Jval = cell2mat(CVresults.val.J(:,1));

%% get ROC curve for ALL predictions:
predVar  = 'ASGRADE_T72';
classThr = 1;
Ytarget = modelData.(predVar)(Jval)>=classThr;
Yact    = activationsVal;
figure

[AUC] = performanceSummaryNeurNet([],Ytarget,Yact,...
    numel(Yact),[],true);
title(sprintf('AUCwhole=%.3g, location=%g',AUC.whole,aa));
%% get CV-AUC
close all
Nsplits = 8;
aa = 1;
% linear index for positions for which predictions have been made:
% set the prediction set to be the set of auscultation areas where atleast
% one position has clean recording:
predVar = 'ASGRADE_T72';
% predVar = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
classThr = 1;
[JvalRedist,Npps] = distributePosCasesEvenly(predictionSet.(predVar)>=classThr,...
                                                               Nsplits);
figures = true;
netName = sprintf('netAR');
pos     = sprintf('pos%g',aa);
AUCfield.(netName).(pos).(predVar).(num2name(classThr)) = zeros(Nsplits,1);

for i=1:Nsplits
    Ytarget = predictionSet.(predVar)(JvalRedist.val{i})>=classThr;
    Yact    = activationsVal(JvalRedist.val{i});
    if figures
        figure
    end
    AUC = performanceSummaryNeurNet([],Ytarget,Yact,...
        numel(Yact),[],figures);
    if figures
        title(sprintf('AUCwhole=%.3g, location=%g',AUC.whole,aa))
    end
    AUCfield.(netName).(pos).(predVar).(num2name(classThr))(i) = AUC.whole;
end

mean(AUCfield.netARG1.pos4.MRGRADE_T72.g1)
computeCImeanEst(AUCfield.(netName).(pos).(predVar).(num2name(classThr)),'3','tDist');


%% test G2 and regression nets against eachother
predVar  = 'ASGRADE_T72';
mean(AUCfield.netMurRegAllPos.(predVar) - AUCfield.netMurG2AllPos.(predVar))
pValue(AUCfield.netMurRegAllPos.(predVar) - AUCfield.netMurG2AllPos.(predVar))

predVar  = 'ARGRADE_T72';
xx = AUCfield.netMurRegPos4.(predVar) - AUCfield.netMurG2Pos4.(predVar);
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
xx = AUCfield.netMurReg.(pos).(predVar).(num2name(classThr)) - ...
     AUCfield.netMurG2.(pos).(predVar).(num2name(classThr)) ;
meanValDiffs(classThr,aa) = mean(xx);
pVals(classThr,aa) = pValue(xx);
end
end


