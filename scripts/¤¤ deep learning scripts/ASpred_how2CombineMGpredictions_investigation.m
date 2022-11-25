% Exploring linear combinations of murmur-grade to improve AS-prediction.
%%
CVresults.val.J
CVresults.val.activations
load CVresults_netMurRegAllPos_valStop_overTrain
%%
[paddedActivMatrix,JatleastOne] = getZeroPaddedActivMatrix(CVresults.val.activations,...
                                                         CVresults.val.J);

weigths = [.5,1,.5,1];
paddedActivMatrix(:,1) = weigths(1)*paddedActivMatrix(:,1);
paddedActivMatrix(:,2) = weigths(2)*paddedActivMatrix(:,2);
paddedActivMatrix(:,3) = weigths(3)*paddedActivMatrix(:,3);
paddedActivMatrix(:,4) = weigths(4)*paddedActivMatrix(:,4);


JnonZero = find(sum(paddedActivMatrix,2)>0);
CVresults.valTot.J   = cell(8,1);
CVresults.trainTot.J = cell(8,1);
for i=1:8
    CVresults.valTot.J{i} = find(CVresults.valTot.I{i});
    CVresults.trainTot.J{i} = find(CVresults.trainTot.I{i});
end
                                                     

AUCvec = zeros(8,1);
Ytruth = cell(8,1);
maxAct = cell(8,1);
for i=1:8
    J = CVresults.valTot.J{i};
    maxAct{i} = sum(paddedActivMatrix(J,:).^3,2);
    Ytruth{i} = HSdata.ASgrade(J)>=1;
    AUCvec(i) = getAUC(Ytruth{i},maxAct{i});
end
mean(AUCvec)

Ytarget = cell2mat(Ytruth);
activ = cell2mat(maxAct);
AUC = getAUC(Ytarget,activ)
figure
performanceSummaryNeurNet([],Ytarget,activ,...
                          [],true);
                      
%%
TargetArray = getTruthCellArray(HSdata,CVresults.val.J,"ASgrade",classThr)

aa = 3;
Ytarget = cell2mat(TargetArray(:,aa));
activ = cell2mat(CVresults.val.activations(:,aa));
AUC = getAUC(Ytarget,activ)

%% Developing model for predicting AVMPG using MG-grade 1-4
[paddedActivMatrix,JatleastOne] = getZeroPaddedActivMatrix(CVresults.val.activations,...
                                                         CVresults.val.J);
X1 = array2table(paddedActivMatrix(JatleastOne,:),'v',{'A','P','T','M'});
Ytrain = array2table(sqrt(HSdata.AVMEANPG_T72(JatleastOne)),'v',{'avmpg'});
noise = array2table(paddedActivMatrix(JatleastOne,:)==0,'v',...
    {'noiseA','noiseP','noiseT','noiseM'});
Yas = array2table(HSdata.ASgrade(JatleastOne)>=1,'v',{'as'});
dataMdl = [X1,X2,X3,noise,Ytrain,Yas];
dataMdl = dataMdl(~isnan(dataMdl.avmpg),:);
Jtrain = intersect(JatleastOne,find(~isnan(dataMdl.avmpg)));
% dataMdl = dataMdl(dataMdl.avmpg>sqrt(5),:);

% glm = fitglm(dataMdl,'avmpg ~ a1^4 + a2 + a3^9 + a4^2','distr','binomial');
% glm = fitglm(dataMdl,'as ~ aor + pul + tri + mit','distr','binomial');
% glm = fitglm(dataMdl,'as ~ x + y + z + w','distr','binomial');
% dataMdl = dataMdl(dataMdl.M>0,:);
% glm = fitglm(dataMdl,'avmpg ~ A');
% glm = fitglm(dataMdl,'avmpg ~ A + A2+ P + T2 + T2*noiseA*noiseP + M*noiseA*noiseP*noiseT');
% glm = fitglm(dataMdl,'avmpg ~ A + P + M + T');
% glm = fitglm(dataMdl,'avmpg ~ A + A:A + P + M + T');
% glm = fitglm(dataMdl,'avmpg ~ A:A + P + M + T');
% glm = fitglm(dataMdl,'avmpg ~ A:A + P + M + T + T:noiseA:noiseP');
% glm = fitglm(dataMdl,'avmpg ~ A:A + P + M + T + T:noiseA:noiseP + M:noiseA:noiseP:noiseT');
formula = 'avmpg ~ A:A + P:P + T + P:P:noiseA + (A:A):noiseP + noiseP:noiseA';
glm = fitglm(dataMdl,formula);
glm.Coefficients;
% glm = fitglm(dataMdl,'avmpg ~ A + A2+ P + M + T2 + T*noiseA*noiseP + M*noiseA*noiseP*noiseT');
% glm = fitglm(dataMdl,'avmpg ~ A + A2+ P + M + T2 + T*noiseA*noiseP + noiseA*noiseP*noiseT:');

Y = dataMdl.avmpg;
Ypred = glm.Fitted.Response;
res = glm.Residuals.Raw;    

close all
% figure
% plot(Ypred,res,'o')

[X,Y,T,AUC] = perfcurve(dataMdl.avmpg>=sqrt(10), glm.predict,true);
plot(X,Y)
display(glm);
title(sprintf('AUC=%g',round(AUC*100,2)))
glm.Coefficients
[glm.Deviance,glm.ModelCriterion.BIC]
%% bootstrapping to one vs all variable predictions
posStr = 'P';
data = dataMdl(dataMdl.(posStr)>0,:);
avmpgThr = 10;
myFunc = @(Iboot) getAUConePosVsAllPos_ASpred(data,Iboot,posStr,avmpgThr);
nIterations = 100;
% compute the CI of correlation between murmur grade and
% aortic-valve-pressure-gradient:
J0 = 1:height(data);
[TestSingleVsAll.(posStr),bootsStat] = bootci(nIterations,...
           {myFunc, J0});
       
TestSingleVsAll.pval = pValue(bootsStat);
%% cross validation to estimate one vs all prediction difference
clear AUC
avmpgThr = 15;
positions = {'A','P','T','M'};
% storage of p-values (which differ for each bootstrap):
% ¤¤ CHOOSE WHIH VARIABLE TO PREDICT; avmpg OR as ¤¤
predVar = "avmpg";
% ¤¤ CHOOSE WHICH POSITION TO COMPARE AGAINST ¤¤
for m=1:4
posStr = positions{m};
% ¤¤ CHOOSE CLASS THRESHOLD ¤¤
% avmpgThr = 10;

NbootExperiments = 5;
for k=1:NbootExperiments
    Nsplits = 8;
    % extract rows where the position has usable audio:
    data0 = dataMdl(dataMdl.(posStr)>0,:);
    n = height(data0);
    if predVar=="avmpg"
        c = cvpartition(n,'KFold',Nsplits);
    else
        c = cvpartition(data0.as>0,'KFold',Nsplits);
    end
    
    for i=1:Nsplits
        % train models on training set:
%         data = data0(c.training(i),:);
        data = data0(c.training(i),:);
        oneVarformula = sprintf('avmpg ~ %s',posStr);
        glm_one = fitglm(data,oneVarformula);
        glm_all = fitglm(data,'avmpg ~ A:A + P + M + T + T:noiseA:noiseP + M:noiseA:noiseP:noiseT');
        
        % predict on validation set:
%         data = data0(c.test(i),:);
        data = data0(c.test(i),:);
        activ_one = glm_one.predict(data);
        activ_all = glm_all.predict(data);
        
        if predVar=="avmpg"
            AUC.(posStr).one(i,1) = getAUC(data.avmpg>=sqrt(avmpgThr), activ_one);
            AUC.(posStr).all(i,1) = getAUC(data.avmpg>=sqrt(avmpgThr), activ_all);
        else
            AUC.(posStr).one(i,1) = getAUC(data.as>=1, activ_one);
            AUC.(posStr).all(i,1) = getAUC(data.as>=1, activ_all);
        end
    end
    p.(predVar).(posStr)(k) = pValue(AUC.(posStr).all-AUC.(posStr).one);
end
p.(predVar).(posStr) = mean(p.(predVar).(posStr));
display(p.(predVar))
end

AUCmat.one.A(avmpgThr-9) = mean(AUC.A.one);
AUCmat.one.P(avmpgThr-9) = mean(AUC.P.one);
AUCmat.one.T(avmpgThr-9) = mean(AUC.T.one);
AUCmat.one.M(avmpgThr-9) = mean(AUC.M.one);

AUCmat.all.A(avmpgThr-9) = mean(AUC.A.all);
AUCmat.all.P(avmpgThr-9) = mean(AUC.P.all);
AUCmat.all.T(avmpgThr-9) = mean(AUC.T.all);
AUCmat.all.M(avmpgThr-9) = mean(AUC.M.all);

%%
close
subplot(2,2,1)
plot(10:15,AUCmat.all.A,'r')
hold on
plot(10:15,AUCmat.one.A,'r','LineStyle','--')
subplot(2,2,2)
plot(AUCmat.all.P,'r')
hold on
plot(AUCmat.one.P,'r')







