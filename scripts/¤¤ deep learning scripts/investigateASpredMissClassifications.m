% Investigating the prediction errors from predicting AS. Run
% estSnSpBalancedWay before running this script.

Jpos = find(Ytrue==1);
Jneg = find(Ytrue==0);
Jpos = JvalJoint(Jpos);
Jneg = JvalJoint(Jneg);

JfalsePos = find(and(predictions==1,Ytrue==0));
JfalsePos = JvalJoint(JfalsePos);
JfalseNeg = find(and(predictions==0,Ytrue==1));
JfalseNeg = JvalJoint(JfalseNeg);

modelData.ASGRADE_T72(JfalsePos)
modelData.ASGRADE_T72(JfalseNeg)

% analyze missclassified cases; false AS predictions.
Jall = 1:height(modelData);
Isub = findInd(Jall,JfalsePos);
meanPGfalsePos = modelData.AVMEANPG_T72(JfalsePos);
meanPGfalseNeg = modelData.AVMEANPG_T72(JfalseNeg);
murGradeFalseNeg = modelData.maxMeanMurGrade(JfalseNeg);
sum(meanPGfalsePos>10)
boxplot(modelData.AVMEANPG_T72,Isub)
mean(meanPGfalsePos,'o')
mean(meanPGfalseNeg,'o')


modelData.LVSTROKEVOL_T72(JfalseNeg).*modelData.HR_ECHO_T72(JfalseNeg)
modelData.LVSTROKEVOL_T72(JfalseNeg)

mean(modelData.BMI_T7(JfalsePos),'o')
mean(modelData.LVSTROKEVOL_T72,'o')








