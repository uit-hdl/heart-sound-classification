% script where I asses annotator interrater agreement.
thr = 2;
IadMax = max([HSdata.MURMUR_1GRADENRAD_T72,HSdata.MURMUR_2GRADENRAD_T72...
    HSdata.MURMUR_3GRADENRAD_T72,HSdata.MURMUR_4GRADENRAD_T72],[],2);
IsaMax = max([HSdata.MURMUR_1GRADENRSA_T72,HSdata.MURMUR_2GRADENRSA_T72...
    HSdata.MURMUR_3GRADENRSA_T72,HSdata.MURMUR_4GRADENRSA_T72],[],2);

Iad = HSdata.MURMUR_1GRADENRAD_T72;
Isa = HSdata.MURMUR_1GRADENRSA_T72;
I3  = HSdata.maxMeanMurGrade;

% [X,Y] = perfcurve(LABELS,SCORES,POSCLASS) computes a ROC curve
[X,Y,T,AUC] = perfcurve(Iad>=thr, Isa, true)
figure
plot(X,Y)

% AUC of predictions across all 4 positions:
% AD predicts SA: 0.9331, 0.9358, 0.9245, 0.9273,
% SA predicts AD: 0.9387, 0.9646, 0.9245, 0.9429

% average AUC across predictions:
% AD predicts SA: 0.9302
% SA predicts AD: 0.9427

% prediction of maximum murmur:
% ADmax predicts max: 0.9892
% SAmax predicts max: 0.9921
%% annotator corellation
for aa=1:4
    HSdata.(sprintf('murG%g_AD',aa)) = HSdata.(sprintf('MURMUR_%gGRADENRAD_T72',aa));
    HSdata.(sprintf('murG%g_SA',aa)) = HSdata.(sprintf('MURMUR_%gGRADENRSA_T72',aa));
end

HSdata.maxMurSA = max([HSdata.murG1_SA,HSdata.murG2_SA,HSdata.murG3_SA,HSdata.murG4_SA],[],2);
HSdata.maxMurAD = max([HSdata.murG1_AD,HSdata.murG2_AD,HSdata.murG3_AD,HSdata.murG4_AD],[],2);

corrEachPos = [corr(HSdata.murG1_AD,HSdata.murG1_SA),...
             corr(HSdata.murG2_AD,HSdata.murG1_SA),...
             corr(HSdata.murG3_AD,HSdata.murG1_SA),...
             corr(HSdata.murG4_AD,HSdata.murG1_SA)]

pooledAnnot_SA = [HSdata.murG1_AD;...
                  HSdata.murG2_AD;...
                  HSdata.murG3_AD;...
                  HSdata.murG4_AD];
pooledAnnot_AD = [HSdata.murG1_SA;...
                  HSdata.murG2_SA;...
                  HSdata.murG3_SA;...
                  HSdata.murG4_SA];  
%%
corr([pooledAnnot_SA,pooledAnnot_AD])

%% cohens kappa
cohensKappa(categorical(pooledAnnot_SA),categorical(pooledAnnot_AD));

kappa = cohensKappa(pooledAnnot_SA>=2,pooledAnnot_AD>=2)

%% prevalence of murmurs>=1

% DIMENSIONS: ANNOTATOR ID; AUSCULTATION POSITION; MURMUR CUTOFF; 

dec2perc(mean(pooledAnnot_SA>=1))
dec2perc(mean(pooledAnnot_AD>=1))

dec2perc(mean(pooledAnnot_SA>=2))
dec2perc(mean(pooledAnnot_AD>=2))

%% bar-plot murmur prevalence across each annotator and each position. 
murPrevMat1 = [[mean(HSdata.murG1_SA>=1),mean(HSdata.murG1_AD>=1)];...
               [mean(HSdata.murG2_SA>=1),mean(HSdata.murG2_AD>=1)];...
               [mean(HSdata.murG3_SA>=1),mean(HSdata.murG3_AD>=1)];...
               [mean(HSdata.murG3_SA>=1),mean(HSdata.murG3_AD>=1)];...
               [mean(HSdata.maxMurSA>=1), mean(HSdata.maxMurAD>=1)]];
murPrevMat2 = [[mean(HSdata.murG1_SA>=2),mean(HSdata.murG1_AD>=2)];...
               [mean(HSdata.murG2_SA>=2),mean(HSdata.murG2_AD>=2)];...
               [mean(HSdata.murG3_SA>=2),mean(HSdata.murG3_AD>=2)];...
               [mean(HSdata.murG3_SA>=2),mean(HSdata.murG3_AD>=2)];...
               [mean(HSdata.maxMurSA>=2), mean(HSdata.maxMurAD>=2)]];

barText = categorical({'aortic','pulmonic','tricuspid','mitral','all'});
barText = reordercats(barText,{'aortic','pulmonic','tricuspid','mitral','all'});

subplot(1,2,1)        
    bar(barText,murPrevMat1*100)
    ylabel 'prevalence (%)'
    ylim([0,25])
    title('murmur-grade $\geq$ 1','Interpreter','latex')

subplot(1,2,2)        
    bar(barText,murPrevMat2*100)
    legend({'annotator: SA','annotator: AD'})
    ylim([0,25])
    title('murmur-grade $\geq$ 2','Interpreter','latex')
    
%% (only non-noisy recordings) bar-plot murmur prevalence across each annotator and each position 
murPrevMat1 = [[condProb(HSdata.murG1_SA>=1,HSdata.noise1==0),condProb(HSdata.murG1_AD>=1,HSdata.noise1==0)];...
               [condProb(HSdata.murG2_SA>=1,HSdata.noise2==0),condProb(HSdata.murG2_AD>=1,HSdata.noise2==0)];...
               [condProb(HSdata.murG3_SA>=1,HSdata.noise3==0),condProb(HSdata.murG3_AD>=1,HSdata.noise3==0)];...
               [condProb(HSdata.murG3_SA>=1,HSdata.noise4==0),condProb(HSdata.murG3_AD>=1,HSdata.noise4==0)];...
               [mean(HSdata.maxMurSA>=1),mean(HSdata.maxMurAD>=1)]];
           
murPrevMat2 = [[mean(HSdata.murG1_SA>=2),mean(HSdata.murG1_AD>=2)];...
               [mean(HSdata.murG2_SA>=2),mean(HSdata.murG2_AD>=2)];...
               [mean(HSdata.murG3_SA>=2),mean(HSdata.murG3_AD>=2)];...
               [mean(HSdata.murG3_SA>=2),mean(HSdata.murG3_AD>=2)];...
               [mean(HSdata.maxMurSA>=2),mean(HSdata.maxMurAD>=2)]];

barText = categorical({'aortic','pulmonic','tricuspid','mitral','combined'});
barText = reordercats(barText,{'aortic','pulmonic','tricuspid','mitral','combined'});

subplot(1,2,1)        
    bar(barText,murPrevMat1*100)
    ylabel 'prevalence (%)'
    ylim([0,25])
    title('murmur-grade $\geq$ 1','Interpreter','latex')

subplot(1,2,2)        
    bar(barText,murPrevMat2*100)
    legend({'annotator: SA','annotator: AD'})
    ylim([0,25])
    title('murmur-grade $\geq$ 2','Interpreter','latex')
    
%% Confidence intervals murmur prevalence
computeCImeanEst(HSdata.murG1_SA>=1,"2")
computeCImeanEst(HSdata.murG2_SA>=1,"2")
computeCImeanEst(HSdata.murG3_SA>=1,"2")
computeCImeanEst(HSdata.murG4_SA>=1,"2")


