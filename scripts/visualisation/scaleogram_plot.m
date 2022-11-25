
k = 0;

%%
J = find(HSdata.murGrade1>=1);

k = k+1;
kk = J(k);
aa = 1;
id = HSdata.id(kk);
[x,fs0] = wav2TS(id,aa);
Nds0 = 20;
x = downsample(x,Nds0);
fs = floor(fs0/Nds0);
close all
% figure('units','normalized','outerposition',[0 0 1 1])
figure
scaleogramPlot(x,fs);
title(sprintf('k=%g',k))