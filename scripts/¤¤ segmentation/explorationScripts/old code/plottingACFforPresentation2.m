% k = k+1;
id = HSdata.UNIKT_LOPENR(k);
aa = 1;
Nds = 20;
NdsAcf = 35;
x = wav2TS(id,aa);
x = downsample(x,Nds);
Fs = floor(44100/Nds);
figures = true;
close all
subplot(211)
runSpringerSegmentationAlgorithm(x, Fs, HMMpar.Bmatrix, HMMpar.piVector, HMMpar.totObsDist, figures)
% modSeg = ngbrSegment(id,Nds,NdsAcf,HMMpar);
% auto corelation function
subplot(212)
[heartRate, systolicTimeInterval] = getHeartRateSchmidt(x, Fs, true);

% clf
% getScaleogram(x,1,true);