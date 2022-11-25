function out = getAUConePosVsAllPos_ASpred(dataMdl,Jboot,posStr,avmpgThr)
% Computes the AUC for prediction of AVMPG >= avmpgThr.
oneVarformula = sprintf('avmpg ~ %s',posStr);
dataMdl = dataMdl(Jboot,:);

% fit multi and single position linear regression models:
glm_one = fitglm(dataMdl,oneVarformula);
glm_all = fitglm(dataMdl,'avmpg ~ A:A + P + M + T + T:noiseA:noiseP + M:noiseA:noiseP:noiseT');

[~,~,~,AUCone] = perfcurve(dataMdl.avmpg>=sqrt(avmpgThr), glm_one.predict,true);
[~,~,~,AUCall] = perfcurve(dataMdl.avmpg>=sqrt(avmpgThr), glm_all.predict,true);
out = AUCall - AUCone;
end