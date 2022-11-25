% in this script I do a first plotting and investigation into the
% publically available dataset I found on:
% https://github.com/yaseen21khan/Classification-of-Heart-Sound-Signal-Using-Multiple-Features-/blob/master/MVP_New_3%EC%A3%BC%EA%B8%B0.rar

%% AS
plotScaleogram = true;
vhd = 'AS';

close all
P = 110;
N = 30;
D = 1;
for i=P:P+N
    S = findSubplotSize(N+1)
    subplot(S(1),S(2),i-P+1)
    ii = i*D;
    if ii<10
        fileName = sprintf('New_%s_00%g.wav',vhd,ii);
    elseif i>=10 && i<100
        fileName = sprintf('New_%s_0%g.wav',vhd,ii);
    else
        fileName = sprintf('New_%s_%g.wav',vhd,ii);
    end
    
    [x,Fs] = audioread(fileName);
    x = downsample(x,4);
    
    if plotScaleogram
        getScaleogram(x,Fs,true);
    else
        plot(x)
    end
end

%% MR
