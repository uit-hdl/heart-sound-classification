function stateVecReconstructed = getExpandedReprOfStates(nodes,stateVals,n)
% reconstructs the state vector uing the compact representation given by
% the nodes (indeces where a new state occurs for the first time),
% state values corresponding to the nodes, and the length of the timeseries
% that corresponds to the state vector. Function is intended to be used in
% conjunction with the function "getCompactReprOfStates".

%%
if iscolumn(nodes)
    nodes = nodes';
end
if iscolumn(stateVals)
    stateVals = stateVals';
end

% test example:
% nodes = [2,4,6,9]
% stateVals = [1,2,3,4]
% nodes = [3 6 7 10]
% stateVals = [4 1 2 3]
% n = 11

preceedingState =@(x) x-1 + (x==1)*4;
initState = preceedingState(stateVals(1));
stateVals = [initState, stateVals];

nodesMod = [1,nodes,n+1];
stateLengths = nodesMod(2:end)-nodesMod(1:end-1);
stateVecReconstructed = zeros(1,n);

for i=1:numel(stateLengths)
    x = repmat(stateVals(i),1,stateLengths(i));
    index = nodesMod(i):(nodesMod(i+1) - 1);
    stateVecReconstructed(index) = x;
%     stateVecReconstructed = [stateVecReconstructed,x];
end

end