function M = getScaleogram(x,Fs,plotIt,crop,ab)
% ** DESCRIPTION **
% computes and displays the scaleogram of the signal x. Fs is the sampling
% frequency of the signal (if downsampled then Fs = Fs0/Nds). 
%% preliminary
if nargin==1
    Fs = 44100;
    plotIt = false;
    crop = 10:60;
    ab = [.85 .03];
elseif nargin==2
    plotIt = false;
    crop = 10:60;
    ab = [.85 .03];
elseif nargin==3
    crop = 10:60;
    ab = [.85 .03];
elseif nargin==4
    ab = [.85 .03];
end
%%



[cfs,~] = cwt(x,'amor',Fs);
% tms = (0:length(x)-1)/Fs; % sampling times
M = abs(cfs);

if crop==false
    crop = 1:height(M);
end

if plotIt
%     a = .85;
%     b = 0.03;
%     ab(1)=1;
%     ab(2)=inf;
    imagesc(pseudoLog(M(crop,:),ab(1),ab(2)))
end

end