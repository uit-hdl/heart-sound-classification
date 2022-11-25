function quickScaleogramPlot(id,aa)
[x,fs] = wav2TS(id,aa);
x = downsample(x,20);
fs = floor(fs/20);
scaleogramPlot(x,fs);
end