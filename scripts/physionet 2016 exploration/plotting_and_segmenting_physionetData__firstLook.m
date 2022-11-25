% read a wav. file from dataset a:
[x,Fs] = audioread('a0002.wav'); %#ok<*ASGLU>


%% plot scaleograms of physionet data:
set = 'e'; 
close all
plotScaleogram = true;
P = 20;
N = 10;
D = 1;
for i=P:P+N
    S = findSubplotSize(N+1)
    subplot(S(1),S(2),i-P+1)
    ii = i*D;
    fileName = getNameOfWav(set,ii);
    
    [x,Fs] = audioread(fileName);
    x = x(1:1*10^4);
    
    if plotScaleogram
        getScaleogram(x,Fs,true);
    else
        plot(x)
    end
    
end

%% segment PCGs, and plot results
[x,Fs] = audioread(getNameOfWav('a',15));

[assignedStates, heartPar, acf, info] = runSpringerSegmentationAlgorithm...
                                           (x, Fs, HMMpar, [], false);
close all
getScaleogram(x,Fs,true);
hold on
plotAssignedStates(assignedStates,[],[],1);

%%
clear N X

% ¤¤¤ SETTINGS ¤¤¤
plotScaleogram = true;
plotMFCC = false;
playSound = true;
set = 'a';

N.ds = 1; % factor by which signal is downsampled
N.cs = 4; % N cycles per input segment
N.os = 2; % N cycles that overlap per pair of segments
N.segPerPCG = 6; % N segments to extract per PCG

% ¤¤¤ SEGMENT ¤¤¤
K = K+1;
[x0,Fs] = audioread(getNameOfWav(set,K));
if playSound
    soundsc(x0(1:5000),Fs)
end

if plotScaleogram
    close all
    getScaleogram(x0(1:15000),Fs,true);
    title(sprintf('%s'),num2class(Ygt.(set).class(K))) %#ok<*CTPCT>
end

[assignedStates, heartPar, acf, info] = runSpringerSegmentationAlgorithm...
                                           (x0, Fs, HMMpar, [], false);
% convert states into a vector showing start of each cycle:
[~,segLines] = states2cycles(assignedStates);

x  = schmidt_spike_removal(x0, Fs);

Nc_available = numel(segLines) - 1; % N cycles available
N.cs_actual  = min(N.cs,Nc_available); % check if the requested number of cycles is available
% the rows of the matrix L contain the start and stop index of each segment
[L,~] = getSegments(segLines,N.os,N.cs_actual,N.segPerPCG); 
Ncomp = [13,200];

X = cell(1,N.segPerPCG);
% cycle through the segments and compute MFCC for each:
for k=1:N.segPerPCG
    % get the segment:
    xk = x(L(k,1):L(k,2));
    % get compact representation of signal in time-frequency domain:
    Mk = getMFCC(xk,floor(2150/N.ds));
    % extract the number of coefficients desired as features:
    Mk = Mk(:,:);
    % normalize MFCC:
    MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
    if ~isempty(Ncomp)
        MkNorm = imresize(MkNorm, Ncomp);
    end
    X{k} = MkNorm;
end 

if plotMFCC
    close all
    imagesc(X{2})
end

% load network to make prediction with:
load networksTrainingSet_valStop.mat

mur = zeros(1,N.segPerPCG);
for i=1:N.segPerPCG
    mur(i) = net.predict(X{i});
end
activation = median(mur)






% t1 = 25e-3;
% t2 = 10e-3;
% t  = 2205^-1; % length of a timestep
% win = hann(floor(t1/t));
% S = stft(xk,"Window",win,"OverlapLength",floor(t2/t),"Centered",false);
% S
% % coeffs = mfcc(S,fs,"LogEnergy","Ignore");
% coeffs = mfcc(S,2205);
% coeffs = coeffs';
                                       
                                       