% in this script I take a look at the predictions of a network to see if
% there are any notable differences between the samples that are correctly
% and incorrectly classified. Assumes a validation set Xval with
% annotations in Yval, and a network net that has been trained on a single
% auscultation location.
%% load net (if nessecary)
load netAA1murG2_v2
%% predict on validation set:
Ypred = classify(net,Xval,'MiniBatchSize',miniBatchSize);

%% see performance curve:
figure()
[AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval,Yval,net,4*Nval,NsegVal,NsegDesired);

%% Change prediction threshold:
% get threshold that achieves peak performance according to performance
% curve:
opt_x_coord = .941176; % get x-coordinate of a good point on performance curve
[~,imin] = min(abs(Y.whole - opt_x_coord))
peakPerformance = [Y.whole(imin),1-X.whole(imin)]
u = T.whole(imin)
%% Make predictions with new prediction threshold:
mur   = categorical(1);
nomur = categorical(0);
Ypred = categorical(double(p.mur>u));
YpredTot = categorical(double(p.murWhole>u*(1 - .02)));

% specificity and sensitivity:
sum(and(YvalTot==YpredTot,YvalTot==nomur))/N.nomur
sum(and(YvalTot==YpredTot,YvalTot==mur))/N.mur
%%
% take most common prediction of data-sample as the prediction
YpredTot = mode(reshape(reshape(Ypred,[Nval,4*NsegVal])',[NsegVal,4*Nval])',2);
YvalTot  = mode(reshape(reshape(Yval,[Nval,4*NsegVal])',[NsegVal,4*Nval])',2);
%%
YpredTot = categorical(zeros(4*Nval,NsegVal));
YvalTot  = categorical(zeros(4*Nval,NsegVal));
for i=1:4
    YpredTot((i-1)*Nval+1 : i*Nval, :) =...
        reshape(Ypred((i-1)*NsegVal*Nval+1 : i*NsegVal*Nval),...
                                                    [Nval,NsegVal]);
    YvalTot((i-1)*Nval+1 : i*Nval, :) =...
    reshape(Yval((i-1)*NsegVal*Nval+1 : i*NsegVal*Nval),...
                                                [Nval,NsegVal]);
end

YpredSeg = YpredTot;
YvalSeg = YvalTot;

YpredTot = mode(YpredTot,2);
YvalTot = mode(YvalTot,2);
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
N.val      = height(dataVal);
N.tot      = height(dataVal)*4;
N.mur      = sum(I.actual.mur);
N.nomur      = sum(I.actual.nomur);
N.falsePos = sum(I.falsePos);
N.falseNeg = sum(I.falseNeg);
N.truePos  = sum(I.truePos); %#ok<*NOPTS>
N.trueNeg  = sum(I.trueNeg)
accuracy = (N.truePos + N.trueNeg)/N.tot



%% investigate missed murmurs (false negatives):
% notation: J contains adresses in the full data set that contains 4*Nval
% number of predictions. JJ contains adresses with respect to the
% validation data set, so numbers between 1 to 212.
J.falseNeg = find(I.falseNeg);
% convert to true index:
JJ.falseNeg = mod(J.falseNeg,Nval);

i_err = 2;
% get auscultation area:
aa = 1 + floor(J.falseNeg(i_err)/Nval);
% find corresponding id-number:
id_err = dataVal.UNIKT_LOPENR(JJ.falseNeg);
% get validation set index:
jj = JJ.falseNeg(i_err);
% get recording using the id-number:
x = wav2TS(id_err(i_err),aa);
x = downsample(x,Nds);

clf
getScaleogram(x,1,true);
murString = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
G = [dataVal.ARGRADE_T72(jj),dataVal.MRGRADE_T72(jj),...
     dataVal.ASGRADE_T72(jj),dataVal.MSGRADE_T72(jj)];
title(sprintf('AR|MR|AS|MS grades = [%g,%g,%g,%g], murGrade=%g, AA=%g',G(1),G(2),G(3),G(4),...
                            dataVal.(murString)(jj),aa))
hold on
scatter(segLinesVal{jj,aa}, 35*ones(1,numel(segLinesVal{jj,aa})),'ro','filled')

playHS(id_err(i_err),aa)
% It seems that most of the ones that were missed were fairly
% understandable. Also, I have not tinkered at all with the murmur
% threshold, and some of these missed murmurs will be captured as we trade
% a little bit of specificity for increased sensitivity. One of the
% recordings was a person with maxed out scores on all 4 diseases, and it
% failed to capture the murmur on 2 out 4 cases where the annotators heard
% murmurs. however, the murmurs were not very strong, so it is
% understandable that these where missed.
%% does it pick up on murmur in any other locations?
% get dataVal indeces for murmur predictions:
JJ.predMur = mod(J.predMur, N.val);
ind_others = find(JJ.predMur==jj)
J.predMur(ind_others)
ceil(J.predMur(ind_others)/Nval)
% in the [4,4,3,3] case, it identifies murmurs in location 1,3 and 4.
% in the case of [0,4,0,0], murmur was identified in locations 1,2 and 3.

% conclusion: in the 2 cases where murmurs were missed, the algorithm
% detected murmurs in the 3 remaining locations for those subjects. The
% algorithm seems to work well.
%% investigate recordings incorrectly labeled as "murmur" (false Positives):
J.falsePos = find(I.falsePos);
% convert to true index:
JJ.falsePos = mod(J.falsePos,Nval);
% find id-number corresponding to recording:
id_err = dataVal.UNIKT_LOPENR(JJ.falsePos);

i_err = 1;
% get auscultation area:
aa = 1 + floor(J.falsePos(i_err)/Nval);
% get ordered index:
jj = JJ.falsePos(i_err);
% get recording using the id-number:
x = wav2TS(id_err(i_err),aa);
x = downsample(x,Nds);

clf
getScaleogram(x,1,true);
murString = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
G = [dataVal.ARGRADE_T72(jj),dataVal.MRGRADE_T72(jj),...
     dataVal.ASGRADE_T72(jj),dataVal.MSGRADE_T72(jj)];
title(sprintf('AR|MR|AS|MS grades = [%g,%g,%g,%g], murGrade=%g, AA=%g',G(1),G(2),G(3),G(4),...
                            dataVal.(murString)(jj),aa))

hold on
scatter(segLinesVal{jj,aa}, 35*ones(1,numel(segLinesVal{jj,aa})),'ro','filled')
% id = 10081212 is a very strange case, it has a sound that sounds highly
% characteristic of severe AS, and yet the echo-labels has this
% characterized as no VHD, only mild MR (grade 1). Very understandably, and
% perhaps desirably, the network classifies this sample as murmur. This
% subject has a mean pressure gradient of 17.9, max velocity of 2.75, and
% aortic valve area of 1.37. The cutoff values for mild AS are >20, >2.5, and
% <1.5, so two criteria are met.
ASechoData = [HSdata.AVMAXPG_T72(HSdata.UNIKT_LOPENR==id),... 
HSdata.AVMEANPG_T72(HSdata.UNIKT_LOPENR==id),... 
HSdata.AVVMAX_T72(HSdata.UNIKT_LOPENR==id),...
-HSdata.AVAVMAX_T72(HSdata.UNIKT_LOPENR==id)]
ASechoData(2:end)>[20,2.5,-1.5]
playHS(id_err(i_err),aa);

%% looping through the false positives:
J.falsePos = find(I.falsePos);
% convert to true index:
JJ.falsePos = mod(J.falsePos,Nval);
% find id-number corresponding to recording:
id_err = dataVal.UNIKT_LOPENR(JJ.falsePos);

G = zeros(N.falsePos,4);
MG = zeros(N.falsePos,1);
plotIt = false;
for i=1:N.falsePos
    % get auscultation area:
    aa = 1 + floor(J.falsePos(i)/Nval);
    % get ordered index:
    jj = JJ.falsePos(i);
    % get recording using the id-number:
    x = wav2TS(id_err(i),aa);
    x = downsample(x,Nds);

    murString = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    MG(i)  = dataVal.(murString)(jj);
    G(i,:) = [dataVal.ARGRADE_T72(jj),dataVal.MRGRADE_T72(jj),...
              dataVal.ASGRADE_T72(jj),dataVal.MSGRADE_T72(jj)];

    
    if plotIt
        clf
        getScaleogram(x,1,true);
        title(sprintf('AR|MR|AS|MS grades = [%g,%g,%g,%g], murGrade=%g, AA=%g',...
                                    G(i,1),G(i,2),G(i,3),G(i,4),...
                                    dataVal.(murString)(jj),aa))
        hold on
        scatter(segLinesVal{jj,aa}, 35*ones(1,numel(segLinesVal{jj,aa})),'ro','filled')
        pause(1)
    end
end

clear errorAnal
% number of people with atleast some degree of AR or MR:
errorAnal.valSet_arORmr = mean(or(dataVal.ARGRADE_T72>0,dataVal.MRGRADE_T72>0));
% number of false positives with atleast some degree of AR or MR:
errorAnal.falsePos_arORmr = mean(or(G(:,1)>0, G(:,2)>0) );

% proportion of people with MR grade 1 or higher:
errorAnal.valSet_mr = mean(dataVal.MRGRADE_T72>0);
% proportion of false positives with MR grade 1 or higher:
errorAnal.falsePos_mr = mean(G(:,2)>0);

% number of people with average murmur greater than 0:
errorAnal.valSet_murGE0 = 1/4*(mean(dataVal.Murmur_1_grade_ref_ny_T72>0)+...
                               mean(dataVal.Murmur_2_grade_ref_ny_T72>0)+...
                               mean(dataVal.Murmur_3_grade_ref_ny_T72>0)+...
                               mean(dataVal.Murmur_4_grade_ref_ny_T72>0))
    
% number of false positives with average murmur greater than 0:
errorAnal.falsePos_murGE0 = mean(MG>0)

errorAnal.falsePos_meanMurGrade = mean(MG);
errorAnal.valSet_meanMurGrade = mean(dataVal.Murmur_1_grade_ref_ny_T72);
%%% old comments %%%
% interestingly, 77% of the false negative cases had MR grade 1 or higher.
% Also, 77% (10 of 13) had average murmur rating greater than 0. In
% contrast, 18.9% of the dataset has average murmur grade higher than 0,
% and 

%%% new comments %%%
% As in the analysis of the location 1 predictions, we see that AR, MR and
% weak murmurs are overrepresented in the false positives. There proportion
% with AR or MR degree atleast one increases by about 20%, the proportion
% with MR degree atleast 1 increases by about 22%, and the proportion with
% weakly audible murmur by atleast one annotator increases from 13% to 59%.

%% looping through the false negatives:
J.falseNeg = find(I.falseNeg);
% convert to true index:
JJ.falseNeg = mod(J.falseNeg,Nval);
% find id-number corresponding to recording:
id_err = dataVal.UNIKT_LOPENR(JJ.falseNeg);

G = zeros(N.falseNeg,4);
MG = zeros(N.falseNeg,1);
plotIt = true;
for i=1:N.falseNeg
    % get auscultation area:
    aa = 1 + floor(J.falseNeg(i)/Nval);
    % get ordered index:
    jj = JJ.falseNeg(i);
    % get recording using the id-number:
    x = wav2TS(id_err(i),aa);
    x = downsample(x,Nds);

    murString = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    MG(i)  = dataVal.(murString)(jj);
    G(i,:) = [dataVal.ARGRADE_T72(jj),dataVal.MRGRADE_T72(jj),...
              dataVal.ASGRADE_T72(jj),dataVal.MSGRADE_T72(jj)];

    
    if plotIt
        clf
        getScaleogram(x,1,true);
        title(sprintf('AR|MR|AS|MS grades = [%g,%g,%g,%g], murGrade=%g, AA=%g',...
                                    G(i,1),G(i,2),G(i,3),G(i,4),...
                                    dataVal.(murString)(jj),aa))
        hold on
        scatter(segLinesVal{jj,aa}, 35*ones(1,numel(segLinesVal{jj,aa})),'ro','filled')
        pause(10)
    end
end

clear errorAnal
% number of people with atleast some degree of AR or MR:
errorAnal.valSet_arORmr = mean(or(dataVal.ARGRADE_T72>0,dataVal.MRGRADE_T72>0));
% number of false positives with atleast some degree of AR or MR:
errorAnal.falseNeg_arORmr = mean(or(G(:,1)>0, G(:,2)>0) );

% proportion of people with MR grade 1 or higher:
errorAnal.valSet_mr = mean(dataVal.MRGRADE_T72>0);
% proportion of false positives with MR grade 1 or higher:
errorAnal.falseNeg_mr = mean(G(:,2)>0);

% number of people with average murmur greater than 0:
errorAnal.valSet_murGE0 = 1/4*(mean(dataVal.Murmur_1_grade_ref_ny_T72>0)+...
                               mean(dataVal.Murmur_2_grade_ref_ny_T72>0)+...
                               mean(dataVal.Murmur_3_grade_ref_ny_T72>0)+...
                               mean(dataVal.Murmur_4_grade_ref_ny_T72>0))

errorAnal.falseNeg_meanMurGrade = mean(MG);
errorAnal.valSet_meanMurGrade = mean(dataVal.Murmur_1_grade_ref_ny_T72)

% Comments: The proportions with AR and/or MR are not significantly
% different from the proportions of the validation set as a whole. Out of
% the nine false negatives, all but 3 were grade 2, with the three being of
% grade 3 or less, and the mean murmur grade being 2.2778. This indicates
% that the missed cases were all borderline or close to borderline cases.

%%% new comments %%%
% there are now only 2 out of 36 murmurs grade 2 or higher that the
% classifier misses, one of grade 2 and one of grade 3. Both these cases
% are of clinical interest, as one has VHD profile [4,4,3,3] (grade 2
% murmur) and the other has VHD profile [0,4,0,0] (grade 3 murmur).
% In both cases, the algorithm detected murmurs in the remaining 3
% locations. This is all the while maintaining a pretty decent specificity
% of 88.4%. It should be noted that the threshold that led to these results
% was somewhat cherrypicked so as to produce a good tradeoff on this
% particular dataset, and we should expect the performance to drop somewhat
% on the test set.

%% MAIN TAKEAWAYS
% AR and particulary MR are overrepresented amongst the false postivies,
% indicating that perhaps the algorithm is capable on picking up on murmurs
% or abnormalities that the annotators did not percieve. Amongst the false
% negatives, they are all borderline or close to borderline cases, and
% diseases are neither under nor overrepresented amongst them. It did
% however miss one very severe case [4,4,3,3], but in that case the
% audibility of the murmur seemed very low, and it was not entirely clear
% why the annotators scored it a 2, as I could neither hear it nor see it.

