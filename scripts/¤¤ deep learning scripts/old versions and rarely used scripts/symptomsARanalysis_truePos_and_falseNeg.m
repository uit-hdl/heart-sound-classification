% first run script estSnSpBalancedWay_valStop to obtain variables for this
% script. Requires TargetMatrix.

% compute ci's for proportions that have symptoms for each correct
% predictions and missed cases:
% clear sym symptomVar
% symptomVar = zeros(height(HSdata),4);
% % symptomVar(1,:) = HSdata.ANGINA_T7;
% % symptomVar(2,:) = HSdata.DYSPNEA_CALMLY_FLAT_T7;
% % symptomVar(3,:) = HSdata.DYSPNEA_FAST_UPHILL_T7; %#ok<*NASGU>
% % symptomVar(4,:) = HSdata.CHEST_PAIN_NORMAL_T7;
% % symptomVar(:,5) = 1;
% symptomVar = sum(symptomVar,2,'o')>0;


sym.ar.pop = computeCImeanEst(symptomVar,"2")*100;
sym.ar.missedPos = computeCImeanEst(symptomVar(JfalseNeg),"2")*100 %#ok<*NOPTS>
sym.ar.caughtPos = computeCImeanEst(symptomVar(JtruePos),"2")*100
sym.ar.missedNeg = computeCImeanEst(symptomVar(JfalseNeg),"2")*100
sym.ar.allTruePos = computeCImeanEst(symptomVar(HSdata.(targetVar)>=classThr),"2")*100
Ivhd = disease2index(targetVarType);
condProb(predMatrix(:,Ivhd),and(symptomVar,targetMatrix(:,Ivhd)))
condProb(predMatrix(:,Ivhd),targetMatrix(:,Ivhd))

close all
plotCIforEachCat([1,2,3],[sym.ar.pop',sym.ar.missedPos',...
                          sym.ar.caughtPos'])
% set(gca,'xtick',1:2,'xticklabel',{'1';'2'});
% boxplot(symptomVar(JfalseNeg))

%% hypotesis test
p1 = sym.ar.missedPos(2)/100;
p2 = sym.ar.caughtPos(2)/100;
n1 = numel(JfalseNeg);
n2 = numel(JtruePos);


% hypothesis: the proportion of missed AR that has symptoms is higher than
% the proportion of the correctly identified cases of AR that has
% symptoms. In other words, the cases that have been caught have a higher
% prevalence of symptoms than the cases that where missed; p2-p1>0
SE = sqrt(p2*(1-p2)/n1 + p1*(1-p1)/n1);
confInt = p2-p1 + [-1 1]*norminv(0.975)*sqrt(p2*(1-p2)/n1 + p1*(1-p1)/n1);
tvalue = p2-p1;
1-normcdf(tvalue/SE)



