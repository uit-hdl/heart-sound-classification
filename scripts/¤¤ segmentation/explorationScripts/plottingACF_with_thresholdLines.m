% Here I plot the ACF along with candidate peaks, main peak, and I plot the
% threshold values used to determine if the tallest peak should be used or
% substituted by a smaller candidate peak for heart rate estimation.

Nds0   = 20;
NdsAcf = 35;
fs0 = 44100;
fs = floor(44100/Nds0);

% get the signal:
k = 6;
freeze = false;
%%
close all
% interesting cases: k = [33, 491, 530]
if not(freeze)
    k=k+1;
end
aa = 1;
id = HSdata.id(k);
x = downsample(wav2TS(id,aa),Nds0);

[hr, sysDur, acFcn, F] = getHeartRateSchmidtOrgl(x, fs, false);
plot((1:numel(acFcn))/fs,acFcn)
hold on
scatter(F.max_index/fs,acFcn(round(F.max_index)),'k*')
scatter(F.true_index/fs,acFcn(round(F.true_index)),'gs')
% scatter(F.index_cand/fs,acFcn(round(F.index_cand)),'gs')
Ix = 1/2*F.max_index/fs*[0.85 1.15];
xline(Ix,'r')
yline(acFcn(round(F.max_index))*0.6,'r')
xlabel 'time-lag (s)'
xline(0.5,'--')
xline(2,'--')
legend({'ACF','maximum','threshold values'})


%%
[assignedStates,heartPar] = runSpringerSegmentationAlgorithm(x, fs, HMMpar);
close all
getScaleogram(x,1,true);
hold on
states2plot = [2,4];
cols        = ['r','g'];
murStr = sprintf(sprintf('Murmur_%g_grade_ref_ny_T72',aa));
plotAssignedStates(assignedStates,states2plot,cols,1,50)
title(sprintf('k=%g, heartRate=%.3g, sys. dur. = %.2g, murGrade=%g',...
    k, heartPar.rate,heartPar.sysDuration, data.(murStr)(k)))




%% junk
acf = getACFschmidt(x,fs);
close all
figure
plot(acf);

