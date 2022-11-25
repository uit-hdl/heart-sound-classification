function [segBorders,kmax] = getSegments(cycleLines,Nos,Ncs,NsegDesired)
% computes endpoints of the segments, each of which contain Ncs cycles and
% overlap with Nos cycles. Extracts NsegDesired segments. If there arent
% enough cycles to meet the desired number, segments are resampled randomly
% untill the number is met. Endpoints are stored as columns in segBorders.
%% Example input
% NsegDesired =6
% Ncs = N.cs
% Nos = N.os
% i=106;
% aa = 3;
% cycleLines = nodes0.seglines{i,aa}
%%
% cycleLines = segLines{i,aa}

% number of cycle lines (lines indicating start of new cycle):
Nl = numel(cycleLines);

if Nl-1<Ncs
    warning('too few cycles')
end

% Number of cycles per segment minus Number of cycles that overlap
x = Ncs - Nos;
% number of segments that can be extracted from the data:
kmax = floor((Nl-Ncs-1)/x + 1);

% compute right endpoints of segments:
rk = Ncs+1+(0:(kmax-1))*x;
% compute left endpoints of segments:
lk = rk-Ncs;
% collect endpoints in a matrix:
segBorders = [lk',rk'];


% If there are not enough segments to meet the desired number, draw
% segments randomly untill condition is met:
NrandSeg = NsegDesired - kmax;

if NrandSeg>0
    ind_randSeg = randi(Nl-Ncs,NrandSeg,1);
    randSeg     = ind_randSeg + [0,Ncs];
    segBorders  = [segBorders;randSeg];
    
else
    segBorders = segBorders((1:NsegDesired),:);
end

% get true index:
segBorders = cycleLines(segBorders);
end