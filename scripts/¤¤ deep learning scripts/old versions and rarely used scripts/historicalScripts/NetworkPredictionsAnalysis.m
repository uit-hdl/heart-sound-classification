% in this script I take a look at the predictions of a network to see if
% there are any notable differences between the samples that are correctly
% and incorrectly classified. Assumes a validation set Xval with
% annotations in Yval, and a network net that has been trained on a single
% auscultation location.
%% load net (if nessecary)
load netAA1murG2_v2
%% predict
Ypred = classify(net,Xval,'MiniBatchSize',miniBatchSize);
%% see performance curve:
figure()
[AUC,X,Y,T] = performanceSummaryNeurNet(Xval,Yval,net,Nval,NsegPerPCG,true);
%% Change prediction threshold

%%
% take most common prediction of data-sample as the prediction
YpredTot = mode(reshape(Ypred,[Nval,NsegVal]),2);
YvalTot  = mode(reshape(Yval,[Nval,NsegVal]),2);
%% whole segments:
dataVal = HSdataVal;
clearvars I N
mur   = categorical(1);
nomur = categorical(0);
I.pred.mur   = (YpredTot==mur)';
I.pred.nomur = (YpredTot==nomur)';
I.actual.mur = (YvalTot==mur)';
I.actual.nomur = (YvalTot==nomur)';


I.truePos  = and(I.pred.mur,I.actual.mur);   % correctly classified murmur
I.falsePos = and(I.pred.mur,I.actual.nomur); % incorrectly classified as murmur
I.trueNeg  = and(I.pred.nomur,I.actual.nomur);
I.falseNeg = and(I.pred.nomur,I.actual.mur); % murmurs that were missed
N.tot     = height(dataVal);
N.falsePos = sum(I.falsePos);
N.truePos = sum(I.truePos); %#ok<*NOPTS>
N.trueNeg = sum(I.trueNeg);
N.falseNeg = sum(I.falseNeg)



%% investigate missed murmurs (false negatives):
i_err  = 3;
id_err = dataVal.UNIKT_LOPENR(I.falseNeg);
k = find(dataVal.UNIKT_LOPENR==id_err(i_err));
x = wav2TS(id_err(i_err),aa);
x = downsample(x,Nds);

clf
getScaleogram(x,1,true);
G = [dataVal.ARGRADE_T72(k),dataVal.MRGRADE_T72(k),...
     dataVal.ASGRADE_T72(k),dataVal.MSGRADE_T72(k)];
title(sprintf('AR|MR|AS|MS grades = [%g,%g,%g,%g], murGrade=%g',G(1),G(2),G(3),G(4),...
                            dataVal.Murmur_1_grade_ref_ny_T72(k)))

hold on
scatter(segLinesVal{k,aa}, 35*ones(1,numel(segLinesVal{k,aa})),'ro','filled')

%% investigate recordings incorrectly labeled as "murmur" (false Positives):
i_err = 13;
% get id numbers of the errors:
id_err = dataVal.UNIKT_LOPENR(I.falsePos);
% get regular index of the error cases:
k = find(dataVal.UNIKT_LOPENR==id_err(i_err));
x = wav2TS(id_err(i_err),aa);
x = downsample(x,Nds);

clf
getScaleogram(x,1,true);
G = [dataVal.ARGRADE_T72(k),dataVal.MRGRADE_T72(k),...
     dataVal.ASGRADE_T72(k),dataVal.MSGRADE_T72(k)];
title(sprintf('AR|MR|AS|MS grades = [%g,%g,%g,%g], murGrade=%g',G(1),G(2),G(3),G(4),...
                            dataVal.Murmur_1_grade_ref_ny_T72(k)))

hold on
scatter(segLinesVal{k,aa}, 35*ones(1,numel(segLinesVal{k,aa})),'ro','filled')
%%
id_err = dataVal.UNIKT_LOPENR(I.falsePos);
G = zeros(N.falsePos,4);
MG = zeros(N.falsePos,1);
plotIt = false;
for i=1:N.falsePos
    x = wav2TS(id_err(i),aa);
    x = downsample(x,Nds);
    k = find(dataVal.UNIKT_LOPENR==id_err(i));
    G(i,:) = [dataVal.ARGRADE_T72(k),dataVal.MRGRADE_T72(k),...
              dataVal.ASGRADE_T72(k),dataVal.MSGRADE_T72(k)];
    MG(i) = dataVal.Murmur_1_grade_ref_ny_T72(k);
    
    if plotIt
        clf
        getScaleogram(x,1,true);
        title(sprintf('AR|MR|AS|MS grades = [%g,%g,%g,%g], murGrade=%g',...
                                    G(i,1),G(i,2),G(i,3),G(i,4),...
                                    MG(i)))
        hold on
        scatter(segLinesVal{k,aa}, 35*ones(1,numel(segLinesVal{k,aa})),'ro','filled')
        pause(1)
    end
end

clear errorAnal
% number of people with atleast some degree of AR or MR:
errorAnal.falsePos.arORmr = mean(or(dataVal.ARGRADE_T72>0,dataVal.MRGRADE_T72>0));
% number of false positives with atleast some degree of AR or MR:
errorAnal.valSet.arORmr = mean(or(G(:,1)>0, G(:,2)>0) );

% number of people with MR grade 1 or higher:
errorAnal.falsePos.mr = mean(dataVal.MRGRADE_T72>0)
% number of false positives with MR grade 1 or higher:
errorAnal.falsePos.mrmean(G(:,2)>0)

% number of people with average murmur greater than 0:
mean(dataVal.Murmur_1_grade_ref_ny_T72>0)
% number of false positives with average murmur greater than 0:
mean(MG>0)


% interestingly, 77% of the false negative cases had MR grade 1 or higher.
% Also, 77% (10 of 13) had average murmur rating greater than 0. In
% contrast, 18.9% of the dataset has average murmur grade higher than 0,
% and 
