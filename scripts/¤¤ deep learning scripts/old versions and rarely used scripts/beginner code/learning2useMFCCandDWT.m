% load audio
data = HSdataTrain;
k = k+1;
id = data.UNIKT_LOPENR(k);
aa = 1;
Nds0 = 20;
x = wav2TS(id,1);
x = downsample(x,Nds0);
fs = 44100/Nds0;
t = fs^-1; % amount of time per timestep
% plot
subplot(311)
getScaleogram(x,1,true);

m = 7;
t1 = 25e-3; % n timesteps in a window
t2 = 10e-3; % number of timesteps overlap between windows
win = hann(floor(t1/t));
S = stft(x,"Window",win,"OverlapLength",floor(t2/t),"Centered",false);
% coeffs = mfcc(S,fs,"LogEnergy","Ignore");
coeffs = mfcc(S,fs);
% nbins = 60;
% coefficientToAnalyze = 4;

% histogram(coeffs(:,coefficientToAnalyze+1),nbins,"Normalization","pdf")
% title(sprintf("Coefficient %d",coefficientToAnalyze))


% plot the mfcc:
subplot(312)
imagesc(flipud(coeffs(:,1:8)'))
axis tight
colormap parula

% get discrete wavelet transform:
[cD cA] = getDWT(x,8,'morl');
subplot(313)
imagesc(abs(cD))


