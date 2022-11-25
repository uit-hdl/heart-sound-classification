% in this script I explore the multiposition model which predicts AS. 

% * develop model

% * Run cross validation to estimate the difference in performance between
% the multiposition model and the single position models.

% * plot results (creates a figure in the article)

%% Developing model for predicting AVMPG using MG-grade 1-4
[paddedActivMatrix,JatleastOne] = getZeroPaddedActivMatrix(CVresults.val.activations,...
                                                           CVresults.val.J);
X = array2table(paddedActivMatrix,'v',{'A','P','T','M'});
Y = array2table(HSdata.avmeanpg,'v',{'avmpg'});
Ysqrt = array2table(sqrt(HSdata.avmeanpg),'v',{'sqrtAvmpg'});
noise = array2table(paddedActivMatrix==0,'v',...
            {'noiseA','noiseP','noiseT','noiseM'});
Yas = array2table(HSdata.ASgrade>=1,'v',{'as'});
dataMdl = [X,noise,Ysqrt,Y,Yas];
Jtrain = intersect(JatleastOne,find(~isnan(dataMdl.avmpg)));
dataMdl = dataMdl(Jtrain,:);

% formula = 'sqrtAvmpg ~ A:A';
% formula = 'avmpg ~ A + A2+ P + T2 + T2*noiseA*noiseP + M*noiseA*noiseP*noiseT');
% formula = 'avmpg ~ A + P + M + T');
% formula = 'avmpg ~ A + A:A + P + M + T');
% formula = 'avmpg ~ A:A + P + M + T');
% formula = 'avmpg ~ A:A + P + M + T + T:noiseA:noiseP');
formula = 'avmpg ~ A:A + P:P + T + (A:A):noiseP + noiseP:noiseA:noiseT:M';
glm = fitglm(dataMdl,formula);
glm.Coefficients.Estimate

Ypred = dataMdl.avmpg>=10;
res   = glm.predict(dataMdl);   

close all
[X,Y,T,AUC] = perfcurve(dataMdl.sqrtAvmpg>=sqrt(10),...
                        glm.predict(dataMdl),true);
plot(X,Y)
display(glm);
title(sprintf('AUC=%g',round(AUC*100,2)))
glm.Coefficients
[glm.Deviance,glm.ModelCriterion.BIC] %#ok<*NOPTS>
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
clear AUC AUCmat pvals
minAvmpg = 7;
maxAvmpg = 30;
for avmpgThr=minAvmpg:maxAvmpg
positions = {'A','P','T','M'};
% storage of p-values (which differ for each bootstrap):
% ¤¤ CHOOSE WHIH VARIABLE TO PREDICT; avmpg OR as ¤¤
predVar = "avmpg";
% ¤¤ CHOOSE WHICH POSITION TO COMPARE AGAINST ¤¤
for m=1:4
posStr = positions{m};
% ¤¤ CHOOSE CLASS THRESHOLD ¤¤
% avmpgThr = 10;

NbootExperiments = 1;
for k=1:NbootExperiments
    Nsplits = 8;
    % extract rows where the position has usable audio:
    data0 = dataMdl(dataMdl.(posStr)>0,:);
%     data0 = dataMdl;
    IposCases = data0.avmpg>=avmpgThr;
    c = cvpartition(IposCases,'KFold',Nsplits);
    
    for i=1:Nsplits
        % train models on training set:
%         data = data0(c.training(i),:);
        data = data0(c.training(i),:);
        oneVarformula = sprintf('avmpg ~ %s',posStr);
        formula = 'avmpg ~ A:A + P:P + T + (A:A):noiseP + noiseP:noiseA:noiseT:M';
        glm_one = fitglm(data,oneVarformula);
        glm_all = fitglm(data,formula);
        
        % predict on validation set:
%         data = data0(c.test(i),:);
        data = data0(c.test(i),:);
        activ_one = glm_one.predict(data);
        activ_all = glm_all.predict(data);
        
        AUC.(posStr).one(i,1) = getAUC(data.sqrtAvmpg>=sqrt(avmpgThr), activ_one);
        AUC.(posStr).all(i,1) = getAUC(data.sqrtAvmpg>=sqrt(avmpgThr), activ_all);
    end
    p(k,m) = pValue(AUC.(posStr).all-AUC.(posStr).one,"twosided")
end
end
p = mean(p,1);
display(p)

pvals.A(avmpgThr-minAvmpg+1) = p(1);
pvals.P(avmpgThr-minAvmpg+1) = p(2);
pvals.T(avmpgThr-minAvmpg+1) = p(3);
pvals.M(avmpgThr-minAvmpg+1) = p(4);

AUCmat.one.A(avmpgThr-minAvmpg+1) = mean(AUC.A.one);
AUCmat.one.P(avmpgThr-minAvmpg+1) = mean(AUC.P.one);
AUCmat.one.T(avmpgThr-minAvmpg+1) = mean(AUC.T.one);
AUCmat.one.M(avmpgThr-minAvmpg+1) = mean(AUC.M.one);

AUCmat.all.A(avmpgThr-minAvmpg+1) = mean(AUC.A.all);
AUCmat.all.P(avmpgThr-minAvmpg+1) = mean(AUC.P.all);
AUCmat.all.T(avmpgThr-minAvmpg+1) = mean(AUC.T.all);
AUCmat.all.M(avmpgThr-minAvmpg+1) = mean(AUC.M.all);
end

%% plot
col = 'k';
figure
subplot(2,2,1)
plot(minAvmpg:maxAvmpg,100*AUCmat.all.A,col)
hold on
plot(minAvmpg:maxAvmpg,100*AUCmat.one.A,col,'LineStyle','--')
plot(minAvmpg+find(pvals.A<0.05)-1,100*AUCmat.all.A(pvals.A<0.05),'k*')
title 'aortic position'
xlabel 'AVPGmean (mm Hg)'

subplot(2,2,2)
plot(minAvmpg:maxAvmpg,100*AUCmat.all.P,col)
hold on
plot(minAvmpg:maxAvmpg,100*AUCmat.one.P,col,'LineStyle','--')
plot(minAvmpg+find(pvals.P<0.05)-1,100*AUCmat.all.P(pvals.P<0.05),'k*')
title 'pulmonic position'
xlabel 'AVPGmean (mm Hg)'

subplot(2,2,3)
plot(minAvmpg:maxAvmpg,100*AUCmat.all.T,col)
hold on
plot(minAvmpg:maxAvmpg,100*AUCmat.one.T,col,'LineStyle','--')
plot(minAvmpg+find(pvals.T<0.05)-1,100*AUCmat.all.T(pvals.T<0.05),'k*')
title 'tricuspid position'
xlabel 'AVPGmean (mm Hg)'

subplot(2,2,4)
plot(minAvmpg:maxAvmpg,100*AUCmat.all.M,col)
hold on
plot(minAvmpg:maxAvmpg,100*AUCmat.one.M,col,'LineStyle','--')
plot(minAvmpg+find(pvals.M<0.05)-1,100*AUCmat.all.M(pvals.M<0.05),'k*')
title 'mitral position'
legend({'AUC multi-position','AUC single-position','sign. difference'})
xlabel 'AVPGmean (mm Hg)'

