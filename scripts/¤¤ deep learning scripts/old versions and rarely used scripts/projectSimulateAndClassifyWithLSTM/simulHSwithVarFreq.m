function simSound = simulHSwithVarFreq(L,T,murAmp,loc,x0)
%% SIMULATE HEART SOUND WITH MURMUR:
% CYCLE LENGTH AND NUMBER OF CYCLES
% T = 4000;          % length of cycle
% L = 20000;         % length of simulated signal
% loc = [.1 .4];
% DETERMINE SHAPE OF SIGNAL
w = [400 400];     % width of each pulse
h = [50,50];       % parameter that determines the width of the pulses
c = loc*T;     % determines the start of each pulse
Aa = [1       1];  % frequency Amplitudes for each pulse
Af = [-.02 -.010];
Hba = .2;           % Base level of max frequency function
Hbf = .018;             % Base level of amplitude function
% CREATE THE FREUENCY AND AMPLITUDE PULSE FUNCTIONS

fmaxVec = pulseFcn2ndDeg(c,w,h,Af,T,L,Hbf,x0);
ampFcn  = pulseFcn2ndDeg(c,w,h,Aa,T,L,Hba,x0);
% plot the max-frequency and amplitude pulse functions
% subplot(211)
% plot(fmaxVec)
% ylim([0,max(fmaxVec)])
% title 'frequency mode function' 
% subplot(212)
% plot(ampFcn)
% ylim([0,max(ampFcn)])
% title 'amplitude function' 

% SIMULATE THE SIGNAL
sigma = 0.000001;
r = 0.92;


% simulate driving noise:
e = sqrt(sigma)*randn(L, 1); 
% we need a loop since the AR(2) coefficients are now time varying
heartSounds = zeros(1,L);

for i=3:L
    % Create the arma:
    Af = fmax2AR2(fmaxVec(i),r);
    heartSounds(i) = -heartSounds(i-1)*Af(2) - heartSounds(i-2)*Af(3) + e(i);
end

% FILTER THROUGH THE AMPLITUDE FUNCTION TO GET THE COMPLETE SIGNAL:
heartSounds = heartSounds.*ampFcn;
% clf
% plot(heartSounds)
% getScaleogram(heartSounds,1,true);
%%
% GENERATE MURMUR
w = 800;
h = 50;
c = .58 *T;
A = [murAmp;.03];
Hb = [.01;.018];

y = simSinglePulse(c,w,h,A,T,L,Hb,x0);
simSound = y + heartSounds;
% clf
% plot(heartSounds)
% clf
% plot(simSound)
% clf
% getScaleogram(y+heartSounds,1,true);
% clf
% getScaleogram(y,1,true);
% title 'signal'
% % ADD LINES WHERE THE NEW PERIODS START
% for k=1:m
% xline(k*T)
% end
%%
% getScaleogram(y,1,true);
% 
% nsamp          = 80;
% window_overlap = 60;
% nfft           = 2^7;
% spectrogram(ymur,hamming(nsamp),window_overlap,[],'yaxis')

end