k=k+1
x0 = wav2TS(HSdataTrain.UNIKT_LOPENR(k),1);
plot(x0)
x = schmidt_spike_removal(x0, fs);
plot(x)
cutoff = 400; % frequenczy in hertz
order = 5;
x = butterworth_low_pass_filter(x,order,cutoff,fs, false);
plot(x)
cutoff = 25; % frequenczy in hertz
order = 5;
x = butterworth_high_pass_filter(x,order,cutoff,fs, false);
plot(x)
x = normalise_signal(x);
plot(x)
Nds = 30;
x = downsample(x,Nds);


clf
subplot(211)
plot(x0)
subplot(212)
plot(x)