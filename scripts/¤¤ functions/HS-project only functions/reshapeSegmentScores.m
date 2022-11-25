function Y_rs = reshapeSegmentScores(Y,NpcgPerLoc,NsegPerPcg,Nloc)
% takes a vector Y containing classification probabilities or
% classifications and reshapes it into a matrix where each row contains the
% values corresponding to PCG signal i of location j. X has length
% Nts*Nloc*NsegPerPcg and is reshaped into a matrix of dimension
% [Npcg*Nloc, NsegPerPcg]. Npcg is the number of PCG signals in the dataset,
% NsegPerPcg is the number of segments per PCG signal, and Nloc is the
% number of locations from which data has been collected.

if nargin==2
    NsegPerPcg = 6;
    Nloc = 1;
elseif nargin==3
    Nloc = 1;
end

if length(Y)< Nloc*NpcgPerLoc*NsegPerPcg
    % operation has already been performed, do nothing.
    Y_rs = Y;
else
    % Yrs = Y reshaped
    NpcgTot = Nloc*NpcgPerLoc;
    Y_rs = zeros(NpcgTot, NsegPerPcg);

    if iscategorical(Y)
        Y_rs = categorical(Y_rs);
    end

    for i=1:Nloc
        Y_rs((i-1)*NpcgPerLoc+1 : i*NpcgPerLoc, :) = ...
            reshape(Y((i-1)*NsegPerPcg*NpcgPerLoc+1 : i*NsegPerPcg*NpcgPerLoc),...
                                                        [NpcgPerLoc,NsegPerPcg]);
    end
end
end