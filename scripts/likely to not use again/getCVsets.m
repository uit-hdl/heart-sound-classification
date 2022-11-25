function getCVsets(N,K,IposCases)
% Splits the data set with linear index set 1:N K times. If Ipos cases is
% supplied, then the positive sampling is done in such a way that the
% positive samples are distributed evenly.
%% example input
N = height(data);
K = 8;
IposCases = data.MSPRESENCE_T72;
%%


Jall = 1:N;
Jpos = find(IposCases);
Jneg = setdiff(Jall,Jpos);
Npos = numel(Jpos);
Nneg = numel(Jneg);
NminPosCasesPerSet = floor(Npos/K);
% how many sets will recieve leftovers after first distribution:
NleftOvers = rem(Npos,K);
% how many sets gets minimum number of positive cases:
NsetsMinPosCases = Npos - NleftOvers;

% start by distributing to the sets who only gets the minimum amount:
JposCasesRemaining = Jpos;
JvalSets.pos = cell(K,1);
for i=1:K
    Jdraw = randsample(JposCasesRemaining,NminPosCasesPerSet);
    JvalSets.pos{i} = Jdraw;
    JposCasesRemaining = setdiff(JposCasesRemaining,Jdraw);
end
% draw which sets gets the left-over cases:
Jlucky = randsample(1:K,NleftOvers);
for i=1:NleftOvers
    JvalSets.pos{Jlucky(i)} = [JvalSets.pos{Jlucky(i)},JposCasesRemaining(i)];
end

JposCasesRemaining = Jpos;
JvalSets.pos = cell(K,1);
for i=1:K
    JvalSets.neg{i} = 
end

end