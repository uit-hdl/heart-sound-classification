for i=1:2
    for aa=1:4
        CVresults.train.J{i,aa} = find(CVresults.train.I{i,aa});
        CVresults.val.J{i,aa} = find(CVresults.val.I{i,aa});
    end
end
%% define container variables
NrowsTot = height(HSdata);


noise_thr = 0.2;
predMatrix = ActMatTrain>=noise_thr;

target = HSdata.avmeanpg(Jval);
ActMatTrain


corr(HSdata.murGrade1(CVresults.val.I{1,1}), CVresults.val.activ{1,1})
%%
aa = 4;
i  = 1;

[m,i_max] = maxk(CVresults.val.activ{i,aa},30);
J_largestNoiseAct = CVresults.val.J{i,aa}(i_max);

close all
figure
for k=1:20
    subplot(4,5,k)
    J = J_largestNoiseAct(k);
    id = HSdata.id(J);
    MG = HSdata.(sprintf('murGrade%g',i))(J);
    noise = HSdata.(sprintf('noise%g',i))(J);
    avmeanpg = HSdata.avmeanpg(HSdata.id==id);
    quickScaleogramPlot(id,aa);
    title(sprintf('[act,AVPG,MG,Noise]=[%g,%g,%g,%g]',...
                    round(m(k)*100),round(avmeanpg,1),MG,noise) )
end
sgtitle(sprintf('20 largest noise activations, aa=%g, i_{cv}=%g',aa,i)) 

%% get largest aortic valve mean pressure gradient from noise pred. samples
% what is the 5 largest AVmeanPG values from the 20 samples with largest
% noise activation?
aa = 1;
i  = 1;

[m,i_max] = maxk(CVresults.val.activ{i,aa},20);
J_largestNoiseAct = CVresults.val.J{i,aa}(i_max);
[m_avmpg,i_max] = maxk(HSdata.avmeanpg(J_largestNoiseAct),5);

close all
figure
for i=1:4
    id = HSdata.id(J_largestNoiseAct(i_max(i)));
    J = find(HSdata.id==id);
    murGrade = HSdata.(sprintf('murGrade%g',aa))(J);
    subplot(2,2,i)
    quickScaleogramPlot(id,aa)
    title(sprintf('avmeanpg=%g, MG=%g',round(HSdata.avmeanpg(J),2),murGrade))
end






