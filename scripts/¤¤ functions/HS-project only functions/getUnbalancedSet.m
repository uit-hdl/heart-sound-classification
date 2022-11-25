function X = getUnbalancedSet(Xsegments,NunBalanced,outputType,NsegPerHS)
% Takes the number of recordings in the unbalanced training or validation
% set, NunBalanced, and returns the original set as it was BEFORE
% resampling to balance the data-set. outputType is either "seg" or
% "whole", and refers to the desired form of the output.

% Xsegments is a cell-array that contains either the predictions or inputs
% corresponding to the segments of the HS-recordings in column-form. It is
% assummed to contain (at the tail of the column) resampled positive cases
% to balance the data-set. 
%% example input
% Xsegments = YvalBal{1};
% NsegPerHS = 6;
% NunBalanced = NumHS.val(1);
% outputType = "whole"
%% preliminary:
if nargin==1
    outputType = "seg";
    NsegPerHS = 6;
    NunBalanced = size(Xsegments,1)/NsegPerHS;
elseif nargin==2
    outputType = "seg";
    NsegPerHS = 6;
elseif nargin==3
    NsegPerHS = 6;
end

if isempty(NunBalanced)
    NunBalanced = size(Xsegments,1)/6;
end
    
%% function code:

Nbalanced = size(Xsegments,1)/NsegPerHS;
X = reshape(Xsegments,[Nbalanced,NsegPerHS]);
X = X(1:NunBalanced,:);
if outputType=="seg"
    X = reshape(X,[NsegPerHS*NunBalanced,1]);
end

end







%%
% Xbal = reshape(X,[N.segPerPCG*NpcgBalanced,1])
% X    = reshape(X(1:Npcg,:),[N.segPerPCG*Npcg,1])
% A = XvalBal{1}; %#ok<*USENS>
% B = Xval{1};
% Nbalanced = size(A,1)/N.segPerPCG;
% 
% A = reshape(A,[Nbalanced,N.segPerPCG]);
% B = reshape(B,[NumHS.val(1),N.segPerPCG]);

% A(1:NumHS.val(1),:)==B



