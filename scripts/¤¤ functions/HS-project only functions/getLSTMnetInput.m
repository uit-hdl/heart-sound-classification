function X = getLSTMnetInput(x0,fs,segLines,N,Ncomp)
% takes audio x0 and a vector of segmentation lines (indicates beginning of
% each cardiac cycle) and returns a cell array X with an MFCC matrix for
% each segment (size Ncomp(1) x Ncomp(2)). N is a structure that specifies
% segment extraction data processing parameters.

% N.ds = factor by which to downsample PCG signal
% N.cs = number of cycles per segment
% N.os = cycle overlap for each consecutive pair of segments
% N.segPerPCG = number of segments to extract from each PCG signal
%% preliminary
if nargin==4
    Ncomp = [13,200];
end

if ~exist('N.ds','var')
    N.ds = 1;
end
if ~exist('N.cs','var')
    N.cs = 4;
end
if ~exist('N.os','var')
    N.os = 2;
end
if ~exist('N.segPerPCG','var')
    N.segPerPCG = 6;
end
%% code

% Remove spikes and downsample:
x  = schmidt_spike_removal(x0, fs);
x  = downsample(x,N.ds);

Nc_available = numel(segLines) - 1;
N.cs_actual  = min(N.cs,Nc_available);
[L,~] = getSegments(segLines,N.os,N.cs_actual,N.segPerPCG);

% ¤¤¤ get MFCC matrix for each segment ¤¤¤
X = cell(1,N.segPerPCG);
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
    X{1,k} = MkNorm;
end

end