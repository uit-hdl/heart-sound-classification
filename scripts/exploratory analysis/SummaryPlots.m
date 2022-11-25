%% plot VHD grade against murmur loudness
pathStr = 'MR';
if pathStr(2)=='S'
    maxGrade = 3;
else
    maxGrade = 4;
end

x = HSdata.(strcat(pathStr,'grade'));
y = HSdata.maxMeanMurGrade;

[meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(x,y)

clf
for i=1:maxGrade
plot(x(Icat{i})+0.05*randn(sum(Icat{i}),1),y(Icat{i}),'.')
hold on
end
plotCIforEachCat(cat,CI)
title(sprintf('max mean murmur grade vs %s severity grade',pathStr))

%% histogram showing how auscultation-murmur-profile looks for different valvular diseases; MURMUR GRADE WEIGHTED
pathStr = 'MS';
if pathStr(2)=='S'
    maxGrade = 3;
else
    maxGrade = 4;
end  
clf
yLim = 0;
for i=1:maxGrade
    gradeStr = sprintf('grade%g',i);
    Iobs  = stats.I.(pathStr).EQ.(gradeStr).any;
    cols = murDataInfo.meanMurGrade.col;
    xx = sum(HSdata{Iobs, cols});
    yLim = max(yLim,max(xx));
    
    subplot(maxGrade,1,i)
        bar(xx)
        ylim([0,max(xx)*1.1])
        title(sprintf('grade = %g',i))
end
for i=1:maxGrade
    subplot(maxGrade,1,i)
    ylim([0,yLim]*1.1)
end
sgtitle(sprintf('%s audibility for each auscultation locations',pathStr)) 

% For AR we see that as grade increases, the relative murmur audibility at
% location 1 decreases. This might be because having AR means you are more
% likely to have calcifications in general, and we are hearing those the
% loudest at location 1 since thats where mumurs are most commonly heard.
%% histogram showing how auscultation-murmur-profile looks for different valvular diseases; YES/NO murmur
pathStr = 'AR';
if pathStr(2)=='S'
    maxGrade = 3;
else
    maxGrade = 4;
end  
clf

yLim = 0;
for i=1:maxGrade
    gradeStr = sprintf('grade%g',i);
    Iobs  = stats.I.(pathStr).EQ.(gradeStr).any;
    cols = murDataInfo.meanMurGrade.col;
    subplot(maxGrade,1,i)
        yLim = max(yLim,max(xx));
        xx = sum(HSdata{Iobs, cols}>0);
        bar(xx)
        ylim([0,max(xx)*1.1])
        title(sprintf('grade = %g',i))
end
for i=1:maxGrade
    subplot(maxGrade,1,i)
    ylim([0,yLim]*1.1)
end
sgtitle(sprintf('%s audibility for each auscultation locations',pathStr)) 


