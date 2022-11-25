function [paddedActivMatrix,JatleastOne] = getZeroPaddedActivMatrix(activationCell,...
                                                                     Jcell,Nrows)

% Returns a 4 column matrix that been padded with zeros in the positions
% where a predictions are missing due to noise. activationCell contains an
% 8 by 4 cell where each cell contains a a vector of indeces and
% activations correspodning to CV-iteration i and auscultation position aa.
% Jcell contains the corresponding indeces.
%% example input
% activationCell = CVresults.val.activations;
% Jcell = CVresults.val.J;
%% preliminary
if nargin==2
    fullSize = false;
else
    fullSize = true;
end
%%


if isempty(Jcell)
    % if empty, assume that predictions have been made for all positions
    sz = size(activationCell);
    Jcell = cell(sz);
    Jall = (1:numel(activationCell{1,1}))';
    for i=1:sz(1)
        for aa=1:sz(2)
            Jcell{i,aa} = Jall;
        end
    end
end
%% code:
Atot = cell(1,4);
Jtot = cell(1,4);
S = zeros(4,2);
for aa=1:4
    Atot{aa} = cell2mat(activationCell(:,aa));
    Jtot{aa} = cell2mat(Jcell(:,aa));
    S(aa,:)  = getSpan(Jtot{aa});
end

% get range of observations for which there are predictions:
S = [1,max(S(:,2))];

if fullSize
    paddedActivMatrix = zeros(Nrows,4);
else
    paddedActivMatrix = zeros(S(2) - (S(1)-1),4);
end

for aa=1:4
    paddedActivMatrix(Jtot{aa},aa) = Atot{aa};
end

% linear indeces for observations with at least one non-noisy recording:
JatleastOne = unionIterated(Jtot);

end