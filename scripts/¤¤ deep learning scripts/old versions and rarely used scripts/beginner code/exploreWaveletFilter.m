
% get audio signal:
data = HSdataTrain;
k = k+1;
id = data.UNIKT_LOPENR(k);
aa = 1;
Nds0 = 20;
x = wav2TS(id,1);
x = downsample(x,Nds0);
% x = normalise_signal(x);

% get the wavelet filters. D stands for decomposition, R stands for
% reconstruction, Lo stands for low pass, and Hi stands for high pass.
[LoD,HiD,LoR,HiR] = wfilters('db6');
%
y = filter(HiR,1,x);

clf
subplot(311)
plot(x)
subplot(312)
plot(y,'g')
subplot(313)
plot(x)
hold on
plot(y,'g')

%% plot spectrograms
clf
subplot(311)
M0 = getScaleogram(x,1,true);
subplot(312)
M1 = getScaleogram(y,1,true);
subplot(313)
Mcorr =  M0-M1;

imagesc(Mcorr-min(Mcorr,[],'all'))

% for the low pass filter, I see no difference at all between the
% scaleograms. I am starting to suspect that when using time-frequency
% domain, the filter does not matter very much. At least I do not think it
% will be of much importance relative to the other features.
