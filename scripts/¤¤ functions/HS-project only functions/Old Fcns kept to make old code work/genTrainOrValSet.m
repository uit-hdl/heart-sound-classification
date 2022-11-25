function [Xbal,Ybal,X,Y] = genTrainOrValSet(dataFrame,Y0,J,N,aa,balance,nodes,Ncomp,fs,posThr)
% this function generates a training or validation set by extracting
% N.segPerPCG segments from each PCG signal, and stores them in a cell
% array of dimension (number of PCG signals) times (N.segPerPCG).

% aa = auscultation location.
% IposCases = index for positive cases.
% Y0 = numerical or logical vector for regression and binary classes
%     respectively.
% J   = index vector that gives which subset of the data to use
% fs = sampling frequency of original signal
% grade = threshold separating positive from negative cases
% data = data frame corresponding to training or validation set
% nodes.seglines = vector containing the nodes separating the cycles
% N.ds = factor by which to downsample PCG signal
% N.cs = number of cycles per segment
% N.os = cycle overlap for each consecutive pair of segments
% N.segPerPCG = number of segments to extract from each PCG signal
% balance = indicator of whether or not to oversample positive cases to
%           balance dataset.

% extract desired samples:
dataFrame = dataFrame(J,:);
Y0        = Y0(J,:);
segLines  = nodes.seglines(J,:);


if nargin==8
    fs = 44100;
    posThr = 0;
elseif nargin==9
    posThr = 0;
end

if isempty(fs)
    fs = 44100;
end

if isempty(N)
    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
end

if islogical(Y0)
    posThr = 1;
end

IposCases = Y0>=posThr;
JposCases = find(IposCases);

Npcg = height(dataFrame);
% Nbalance is number of additional resamplings that is done to balance the
% dataset, so that Nmur==NnoMur:
if balance
    Nbalance = sum(~IposCases) - sum(IposCases);
else
    Nbalance = 0;
end
% size of inflated data set (added samples to account for few positive samples):
NpcgBalanced = Npcg + Nbalance;

X = cell(NpcgBalanced,N.segPerPCG);
Y = zeros(NpcgBalanced,N.segPerPCG);

for i=1:NpcgBalanced
%     i  %#ok<*NOPRT>
    if i<=Npcg
        ii = i;
        id = dataFrame.UNIKT_LOPENR(ii);
    else
        % randomly pick a sample of category 1 (done to balance data set)
        ii = randsample(JposCases,1);
        id = dataFrame.UNIKT_LOPENR(ii);
    end

    % get time series:
    x0 = wav2TS(id,aa);
    x  = schmidt_spike_removal(x0, fs);
    x  = downsample(x,N.ds);
    
    Nc_available = numel(segLines{ii,aa}) - 1;
    N.cs_actual  = min(N.cs,Nc_available);
    [L,~] = getSegments(segLines{ii,aa},N.os,N.cs_actual,N.segPerPCG);
    
    % cycle through the segments:
    for k=1:N.segPerPCG
        % get the segment:
        xk = x(L(k,1):L(k,2));
        % get compact representation of signal in time-frequency domain:
        Mk = getMFCC(xk,floor(fs/N.ds));
        % extract the number of coefficients desired as features:
        Mk = Mk(:,:);
        % normalize MFCC:
        MkNorm = (Mk-mean(Mk,'all'))/std(Mk,0,'all'); 
        if ~isempty(Ncomp)
            MkNorm = imresize(MkNorm, Ncomp);
        end
        X{i,k} = MkNorm;
        Y(i,k) = Y0(ii);
    end 
end

if islogical(Y0)
    Y = categorical(Y);
end

% reshape to a "vector" form:
Xbal = reshape(X,[N.segPerPCG*NpcgBalanced,1]);
Ybal = reshape(Y,[N.segPerPCG*NpcgBalanced,1]);
X = reshape(X(1:Npcg,:),[N.segPerPCG*Npcg,1]);
Y = reshape(Y(1:Npcg,:),[N.segPerPCG*Npcg,1]);

end