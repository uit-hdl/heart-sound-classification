function z1 = simulHeartSound(murAmp,T,L,loc)

AA = [1 0.95];
CC = 1 ;
arma1 = idpoly( AA, [] , CC );   % constructs an arma process

sigma2 = 0.001;
% L = 5000;
L = L+50;
e1 = sqrt(sigma2)*randn(L, 1) ; % driving noise
y1 = filter(arma1.c , arma1.a , e1) ;% yt = C(z)/A(z)*et
y1 = y1(51:end); % remove samples so that process reaches stationary distribution

% define amplitude function
% T = 900;
% loc = [.1,.3,.45]
t1   = loc(1)*T;
t2   = loc(2)*T;
tMur = loc(3)*T;
A1 = 20; 
A2 = 20;
% murAmp = 15
A3 = murAmp; % murmur amplitude
sigmaAmp = 4*3; % determines the width of the amplitude peaks
sigmaMur = 15*3;
Amp1 = @(x) 0.05+ A1*normpdf(mod(x,T),t1,sigmaAmp)...
                + A2*normpdf(mod(x,T),t2,sigmaAmp)...
                + A3*normpdf(mod(x,T),tMur,sigmaMur);

xplot = 1:L-50;

z1 = Amp1(xplot).*y1';
% plot(z1)
end