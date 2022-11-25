function x1 = preProcessSignal(x0,fs,Nds,passWidth,order)
% preprocess the signal. Removes spikes, filters away high and low
% frequencies using a butterworth filter of order "order", normalizes the
% signal and finally donwsamples it by a factor of Nds.
x = schmidt_spike_removal(x0, fs);
x = butterworth_low_pass_filter(x,order,passWidth(2),fs, false);
x = butterworth_high_pass_filter(x,order,passWidth(1),fs, false);
x = normalise_signal(x);
x1 = downsample(x,Nds);
end