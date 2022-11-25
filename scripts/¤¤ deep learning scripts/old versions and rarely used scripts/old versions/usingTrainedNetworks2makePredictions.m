%% load network:
load netAA1234murG2_v1.mat
load netAA1murG2_v1.mat
load segLinesVal
%% generate segments validation set to test network on
data = HSdataVal;
Nval = height(data);
N.ds = 20;
N.cs = 4;
N.os = 2;
N.segPerPCG = 6;

clf
for i=1:2
    if i==1
        load netAA1murG2_v1.mat
    else
        load netAA1234murG2_v2.mat
    end
    
for aa=1:4
% IposCases = HSdataVal.ARGRADE_T72>=3;
murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = HSdataVal.(murSt)>=2;
balance = false;
Ncomp = [14+16*(i==1),200];
[Xval,Yval] = genTrainOrValSet(data,N,aa,IposCases,balance,segLinesVal,Ncomp);
% predict and show performance curve
miniBatchSize = 2^5;

Ypred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
YpredTot = mode(reshapeSegmentScores(Ypred,Nval,NsegPerPCG,1),2);
YvalTot  = mode(reshapeSegmentScores(Yval, Nval,NsegPerPCG,1),2);
sum(YpredTot(1:Nval)==YvalTot(1:Nval))/(1*Nval)%#ok<*NOPTS>
% sum(Ypred==Yval)/(Nval*6)

[AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval,Yval,net,Nval,N.segPerPCG,false,1);

subplot(4,2,(i-1)*4 + aa)

    plot(Y.seg,1-X.seg)
    hold on
    plot(Y.whole, 1-X.whole)
    plot([0,1],[1,0],'g')
    title(sprintf('AUCseg=%.3g, AUCwhole=%.3g',AUC.seg,AUC.whole))
    legend({'segments','whole PCG-signal'})
    xlabel 'sensitivity'
    ylabel 'specificity'

    pause(.1)
end
end
%%
% take most common prediction of data-sample as the prediction
YpredTot = mode(reshape(Ypred,[Nval,NsegVal]),2);
YvalTot  = mode(reshape(Yval,[Nval,NsegVal]),2);
Nval = height(data);

mur   = categorical(1);
noMur = categorical(0);
accuracy = sum(YpredTot==YvalTot)/Nval

sensTot = sum(and(YpredTot == YvalTot,YpredTot==mur))/sum(YvalTot==mur)
specTot = sum(and(YpredTot==YvalTot,YpredTot==noMur))/sum(YvalTot==noMur)
