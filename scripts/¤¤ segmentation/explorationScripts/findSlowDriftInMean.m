plotAutoCorrelationFcn
%%
% smoothen with butterworth filter
nDs = 3;    
x = downsample(acFcn{3},nDs);
Cbutter = butterworth_low_pass_filter(x,2,1,fs, false);

clf
subplot(211)
plot(x)
xline([1/2,2]*fs/3)
legend({'original signal','smoothed average'})
hold on
plot(Cbutter)
subplot(212)
plot(x-Cbutter)
xline([1/2,2]*fs/3)
legend({'signal minus timevarying mean'})


%% smoothen with gaussian
nDs = 35;
x = downsample(acFcn{3},nDs);
C = smoothdata(x,'gaussian',3200/nDs);
clf
subplot(211)
[pks,locs,w,p] = findpeaks(x(1:110));
plot(x)
xline([1/2,2]*fs/nDs)
hold on
scatter(locs,pks)
legend({'original signal','smoothed average'})
hold on
plot(C)
subplot(212)
plot(x-C)
xline([1/2,2]*fs/nDs)
legend({'signal minus timevarying mean'})

% this method of smoothing seems to work the best.

%% smoothen with homo
homomorphicEnvelope = Homomorphic_Envelope_with_Hilbert(x, Fs);
x = acFcn{3};
C = Homomorphic_Envelope_with_Hilbert(x, Fs);
clf
subplot(211)
plot(x)
legend({'original signal','smoothed average'})
hold on
plot(C)
subplot(212)
plot(x-C)
legend({'signal minus timevarying mean'})

% did not work very well. It compute the upper envelope, so it doesnt work
% well to capture the mean.