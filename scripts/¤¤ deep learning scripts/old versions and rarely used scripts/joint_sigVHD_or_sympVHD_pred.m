% *** script description ***

% in this scipt, I look at the joint prediction of VHD. I consider cases
% prediction of disease as determined by thresholds, as well as prediction
% of symptomatic or assymptomatic disease.

%% prediction of significant and symptomatic VHD:
% to run this script you need: predMatrix and sympAndSickMat

Jpred = unionIterated(CVresults.valTot.J);

% joint predictions
predAtleastOne = sum(predMatrix,2)>0;
sigSickAtleastOne = sum(sympAndSickMat,2)>0;
sigSickAndNoSymp = sum(noSympAndSickMat,2)>0;

Ypred = predAtleastOne;
Ytarget = sigSickAtleastOne;

sympVHD.all.sn = condProb(predAtleastOne(Jpred),sigSickAtleastOne(Jpred));
sympVHD.all.sp = condProb(~predAtleastOne(Jpred),~sigSickAtleastOne(Jpred));

% *** individual predictions ***
% how many detected:
sympVHD.AR.Ncaught = sum(and(predAtleastOne(Jpred),sympAndSickMat(Jpred,1)));
sympVHD.MR.Ncaught = sum(and(predAtleastOne(Jpred),sympAndSickMat(Jpred,2)));
sympVHD.AS.Ncaught = sum(and(predAtleastOne(Jpred),sympAndSickMat(Jpred,3)));
sympVHD.MS.Ncaught = sum(and(predAtleastOne(Jpred),sympAndSickMat(Jpred,4)));

% how many not detected:
sympVHD.AR.Nmissed = sum(and(~predAtleastOne(Jpred),sympAndSickMat(Jpred,1)));
sympVHD.MR.Nmissed = sum(and(~predAtleastOne(Jpred),sympAndSickMat(Jpred,2)));
sympVHD.AS.Nmissed = sum(and(~predAtleastOne(Jpred),sympAndSickMat(Jpred,3)));
sympVHD.MS.Nmissed = sum(and(~predAtleastOne(Jpred),sympAndSickMat(Jpred,4)));

sympVHD.all

%% prediction of significant VHD:
YpredSigVHD = sum(predMatrix,2)>0;
YtargetSigVHD = sum(targetMatrix,2)>0;
sigVHD.all.sn = condProb(YpredSigVHD(Jpred),YtargetSigVHD(Jpred));
sigVHD.all.sp = condProb(~YpredSigVHD(Jpred),~YtargetSigVHD(Jpred));
sigVHD.all

sigVHD.AR.sn = condProb(YpredSigVHD(Jpred),targetMatrix(Jpred,1));
sigVHD.MR.sn = condProb(YpredSigVHD(Jpred),targetMatrix(Jpred,2));
sigVHD.AS.sn = condProb(YpredSigVHD(Jpred),targetMatrix(Jpred,3));
sigVHD.MS.sn = condProb(YpredSigVHD(Jpred),targetMatrix(Jpred,4));

sigVHD.AR
sigVHD.MR
%% difference in symptomaticness; detected vs missed cases
Jpred = unionIterated(CVresults.valTot.J);
sym.(targetVarType).pop = computeCImeanEst(symptomVec,"2")*100;
sym.(targetVarType).missedPos = computeCImeanEst(symptomVec(JfalseNeg),"2")*100 %#ok<*NOPTS>
sym.(targetVarType).caughtPos = computeCImeanEst(symptomVec(JtruePos),"2")*100
sym.(targetVarType).missedNeg = computeCImeanEst(symptomVec(JfalseNeg),"2")*100
sym.(targetVarType).allTruePos = computeCImeanEst(symptomVec(HSdata.(targetVar)>=classThr),"2")*100
Ivhd = disease2index(targetVarType);
IconfirmedVHDandSymp = sympAndSickMat(:,Ivhd);
Ypred = predMatrix(:,Ivhd);
condProb(Ypred(Jpred),IconfirmedVHDandSymp(Jpred))

condProb(~Ypred(Jpred),~IconfirmedVHDandSymp(Jpred))

sym.(targetVarType)

close all
plotCIforEachCat([1,2,3],[sym.(targetVarType).pop',...
                          sym.(targetVarType).missedPos',...
                          sym.(targetVarType).caughtPos'])
                      
%% hypotesis test
p1 = sym.(targetVarType).missedPos(2)/100;
p2 = sym.(targetVarType).caughtPos(2)/100;
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
