%% SIMULATE HEART SOUND WITH MURMUR:
% CYCLE LENGTH AND NUMBER OF CYCLES
m = 10;            % number of cycles to plot
T = 2500;          % length of cycle

% DETERMINE SHAPE OF SIGNAL
w = [300 300 850];   % width of each pulse
h = [25,25,50];      % parameter that determines the width of the pulses
c = [0 700 900];     % determines the start of each pulse
Af = [.08 .12 10]/100;   % frequency Amplitudes for each pulse
Am = 0.5*[90 90 30];
L = 0.02;            % Base level of max frequency function

% CREATE THE FREUENCY AND AMPLITUDE PULSE FUNCTIONS
fmaxVec = pulseFcn(c,w,h,Af,T,L,m);
ampFcn  = pulseFcn2ndDeg(c,w,h,Am,T,1,m);

% plot the max-frequency and amplitude pulse functions
subplot(211)
plot(fmaxVec)
title 'frequency mode function' 
subplot(212)
plot(ampFcn)
title 'amplitude function' 

% SIMULATE THE SIGNAL
sigma = 0.0001;
r = 0.9;

% simulate driving noise:
e = sqrt(sigma)*randn(m*T, 1); 

% we need a loop since the AR(2) coefficients are now time varying
ymur = zeros(1,m*T);

for i=3:T*m
    % Create the arma:
    Af = fmax2AR2(fmaxVec(i),r);
    ymur(i) = -ymur(i-1)*Af(2) - ymur(i-2)*Af(3) + e(i);
end

% FILTER THROUGH THE AMPLITUDE FUNCTION TO GET THE COMPLETE SIGNAL:
y = ymur.*ampFcn;

clf
plot(y)
title 'signal'
% ADD LINES WHERE THE NEW PERIODS START
for k=1:m
xline(k*T)
end
%%
getScaleogram(y,1,true);

signal = y;
[cA,cD] = dwt(signal,'sym4');

lev   = 5;
wname = 'db1'; 
nbcol = 32; 
[c,l] = wavedec(signal,lev,wname);
len = length(signal);
cfd = zeros(lev,len);
for k = 1:lev
    d = detcoef(c,l,k);
    d = d(:)';
    d = d(ones(1,2^k),:);
    cfd(k,:) = wkeep(d(:)',len);
end
cfd =  cfd(:);
I = find(abs(cfd)<sqrt(eps));
cfd(I) = zeros(size(I));
cfd    = reshape(cfd,lev,len);
cfd = wcodemat(cfd,nbcol,'row');

h211 = subplot(2,1,1);
h211.XTick = [];
plot(signal,'r'); 
title('Analyzed signal.');
ax = gca;
ax.XLim = [1 length(signal)];
subplot(2,1,2);
colormap parula;
image(cfd);
tics = 1:lev; 
labs = int2str(tics');
ax = gca;
ax.YTickLabelMode = 'manual';
ax.YDir = 'normal';
ax.Box = 'On';
ax.YTick = tics;
ax.YTickLabel = labs;
title('Discrete Transform, absolute coefficients.');
ylabel('Level');