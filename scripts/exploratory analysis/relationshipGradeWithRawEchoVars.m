
%% MR
x = HSdata.MRgrade;
% x = HSdataRed.maxMeanMurGrade;
y = HSdata.MVREGAREA_T72;
y = g(y,1,1/100);
[meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(x,y);
clf
plot(x + 0.05*randn(numel(x),1), y,'*')
hold on
plotCIforEachCat(cat,CI,'g')
% High correlation:
% MRMAXPG_T72   - Mitral regurgitation maximum pressure gradient
% MRVMAX_T72    - Mitral regurgitation max velocity
% MRVMAX_T72    - Mitral regurgitation flow maximum velocity 
% MRVTI_T72     - Mitral regurgitation flow velocity-time integral
% MVREGAREA_T72 - Mitral Valve Regurgitation Area (strongest by far)
%% MR
x = HSdataRed.ARGRADE_T72;
% x = HSdataRed.maxMeanMurGrade;
y = HSdataRed.MVREGAREA_T72;
y = g(y,1,1/100);
[meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(x,y);
clf
plot(x + 0.05*randn(numel(x),1), y,'*')
hold on
plotCIforEachCat(cat,CI,'g')
% High correlation:
% MRMAXPG_T72   - Mitral regurgitation maximum pressure gradient
% MRVMAX_T72    - Mitral regurgitation max velocity
% MRVMAX_T72    - Mitral regurgitation flow maximum velocity 
% MRVTI_T72     - Mitral regurgitation flow velocity-time integral
% MVREGAREA_T72 - Mitral Valve Regurgitation Area (strongest by far)
%% AS
te.x = HSdataRed.ASGRADE_T72;
te.names = varNamesCell(233:241);
% te.y{1} = HSdataRed.AVMAXPG_T72;
% te.y{2} = HSdataRed.AVAIVMAX_T72;
% te.y{3} = HSdataRed.AVAVMAX_T72;
% te.y{4} = HSdataRed.AVAVTI_T72;
% te.y{5} = HSdataRed.AVCUSP_T72;
% te.y{6} = HSdataRed.AVMEANPG_T72;
% te.y{7} = HSdataRed.AVVMAX_T72;
% te.y{8} = HSdataRed.AVVMEAN_T72;

te.n = height(HSdataRed);

clf
for i=1:8
    subplot(4,2,i)
    te.name = te.names{i};
    te.y = HSdataRed.(te.name);
    [meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(te.x,te.y);

    scatter(te.x + 0.01*randn(te.n,1),te.y)
    hold on
    scatter(cat,meanVals,[],'k')
    
    for j=1:length(cat)
        line([cat(j),cat(j)],[CI(1,j),CI(2,j)],'LineWidth',2,'Color','k')
    end

    title(sprintf('%s',te.names{i}))

end
clearvars te
%% AR
te.x = HSdataRed.ARGRADE_T72;
kk=196+12*2;
te.names = varNamesCell(kk:kk+11);

% te.y{1} = HSdataRed.AVMAXPG_T72;
% te.y{2} = HSdataRed.AVAIVMAX_T72;
% te.y{3} = HSdataRed.AVAVMAX_T72;
% te.y{4} = HSdataRed.AVAVTI_T72;
% te.y{5} = HSdataRed.AVCUSP_T72;
% te.y{6} = HSdataRed.AVMEANPG_T72;
% te.y{7} = HSdataRed.AVVMAX_T72;
% te.y{8} = HSdataRed.AVVMEAN_T72;

te.n = height(HSdataRed);

clf
for i=1:size(te.names,1)
    subplot(4,3,i)
    te.name = te.names{i};
    te.y = HSdataRed.(te.name);
    [meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(te.x,te.y);

    scatter(te.x + 0.01*randn(te.n,1),te.y)
    hold on
    scatter(cat,meanVals,[],'k')
    
    for j=1:length(cat)
        line([cat(j),cat(j)],[CI(1,j),CI(2,j)],'LineWidth',2,'Color','k')
    end

    title(sprintf('%s',te.names{i}))

end
clearvars te

%% AS
x = HSdataRed.ASGRADE_T72;

y = HSdataRed.AVMEANPG_T72;
y = HSdataRed.AVAVMAX_T72;
y = HSdataRed.AVMEANPG_T72;
y = HSdataRed.AVMEANPG_T72./HSdataRed.LVSTRVOLBIP_T72	%./HSdataRed.AVAVMAX_T72;



[meanVals, Icat, yCat, CI, cat] = findMeanForEachCat(x,y);
boxplot(y,x)

figure
plot(x + 0.05*randn(numel(x),1), y,'*')
hold on
plotCIforEachCat(cat,CI,'g')
%%
x = HSdataRed.MSGRADE_T72;

y = HSdataRed.AVMEANPG_T72;
y = HSdataRed.AVAVMAX_T72;
y = HSdataRed.AVMEANPG_T72;
y = HSdataRed.LVOTMEANPG_T72;



[meanVals, Icat, yCat, CI, cat] = findMeanForEachCat(x,y);
clf
plot(x + 0.05*randn(numel(x),1), y,'*')
hold on
plotCIforEachCat(cat,CI,'g')


