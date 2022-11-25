% in this script I want to look for ways to determine whether or not a
% segment should be used for training or if it should be classified as
% "unsegmentable". Including such a category might improve the algorithm,
% as there will be fewer bad samples to confuse it.
%% preliminary
load segLines.mat
load cycles.mat
data = HSdataTrain;
%%

Nds0 = 20;
Fs   = 44100;
k = k+1;
aa = 1;
id = data.UNIKT_LOPENR(k);
x  = wav2TS(id,1);
x = schmidt_spike_removal(x, Fs);   
x = downsample(x,Nds0);
% fpass = 4;
% x = lowpass(x,fpass,Fs/Nds0);
fs = floor(44100/Nds0);
[assignedStates,heartPar,acf,acfInfo] = runSpringerSegmentationAlgorithm(x, fs, HMMpar, []);
[assignedStatesOrgl,heartParOrgl] = runSpringerSegmentationAlgorithmOrgl(x, fs,HMMpar);


clf
subplot(211)
M=getScaleogram(x,1,true);
hold on
states2plot = [2,4];
cols        = ['r','g'];    
murStr = sprintf(sprintf('Murmur_%g_grade_ref_ny_T72',aa));
plotAssignedStates(assignedStates,states2plot,cols,1,50)
plotAssignedStates(assignedStates,states2plot,cols,1,5)
confidence = getSegmentabitityScore(acf,acfInfo);
title(sprintf('k=%g, murGrade=%g, confidence=%.2g, peakProminence=%.2g, peakRatio=%.2g, HRnew/HRold=%.2g',...
              k,data.(murStr)(k), confidence, acfInfo.pp, acfInfo.peakRatio,heartPar.rate/heartParOrgl.rate))
subplot(212)
M=getScaleogram(x,1,true);
title(sprintf('noiseLvl=%.2g',data.MURMUR_1NOISE_REF_T72(k)))
hold on
plotAssignedStates(assignedStatesOrgl,states2plot,cols,1,50)
plotAssignedStates(assignedStatesOrgl,states2plot,cols,1,5)

%%

aa = 1;
id = data.UNIKT_LOPENR(k);
x = wav2TS(id,1);
x = downsample(x,Nds0);
fpass = 80;
xds = lowpass(x,fpass,Fs/Nds0);
subplot(211)
plot(x)
subplot(212)
plot(xds)

lowpass(x,fpass,Fs)
%% spike removal
aa = 1;
id = data.UNIKT_LOPENR(k);
x = wav2TS(id,aa);
[xDespiked] = schmidt_spike_removal(x, Fs);

clf
subplot(211)
plot(x)
subplot(212)
plot(xDespiked)


