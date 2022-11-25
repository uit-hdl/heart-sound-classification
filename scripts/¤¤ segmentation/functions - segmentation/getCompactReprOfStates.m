function [nodes,nodeStates,segLines] = getCompactReprOfStates(stateVec)
% get compact representation of state vector. nodes indicates positions
% where the chain jumps to a new state, i.e. the jump points. nodeStates is
% a vector with state values corresponding to the jump points in nodes.
% segLines indicates positions where a new cycle begins (end of diastole).

nodes = find(stateVec(2:end)-stateVec(1:end-1)) + 1;
nodeStates = stateVec(nodes);
segLines = nodes(nodeStates==1);
end