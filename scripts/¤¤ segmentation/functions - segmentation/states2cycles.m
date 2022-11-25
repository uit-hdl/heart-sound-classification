function [cycleIndex,segLines] = states2cycles(assignedStates)
% converts a vector of states (ex: [2,2,2,3,3,3,...]) to a vector of
% elements indicating heart cycle index. segLines is an index vector
% indicating where each new cycle starts. The start of a cycle is defined
% as the end of diastole.

n = length(assignedStates);
cycleIndex = 1:n;
segLines = [];
cycle = 1;

for i=2:n
    delta = assignedStates(i)-assignedStates(i-1);
    if delta<0
       segLines(end+1) = i; %#ok<AGROW>
       cycle = cycle + 1;
    end
    cycleIndex(i) = cycle;
end

end