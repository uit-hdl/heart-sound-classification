function [J,Npps] = distributePosCasesEvenly(I,Nsplits)
% Function is intended to be used in cross validation to ensure that there
% is equal number of positive cases in each validation set. I is a logical
% index vector that indicates lindear index of positive cases. Returns a
% cell array of size 1 by Nsets that contains linear indeces for each
% validation set after ensuring equal numbers of positive cases in each.
%% example input
% Nsplits = 8;
% I = rand(100,1)>0.7;
%%
N = numel(I);
Jall = 1:N;
Jpos = find(I);
Jneg = setdiff(Jall,Jpos)';
Npos = numel(Jpos);
Nneg = numel(Jneg);
% number of positive examples per set:
Npps = floor(Npos/Nsplits);
Nnps = floor(Nneg/Nsplits);
remainderPos = rem(Npos,Nsplits);
remainderNeg = rem(Nneg,Nsplits);
JvalRedist = cell(Nsplits,2);

xx= ones(1,Nsplits)*Npps;
xx(1:remainderPos) = xx(1:remainderPos) + 1;
xx = cumsum(xx);
xx = [0 xx];
for i=1:Nsplits
    JvalRedist{i,1} = Jpos( xx(i)+1:xx(i+1) );
end

xx = ones(1,Nsplits)*Nnps;
xx(1:remainderNeg) = xx(1:remainderNeg) + 1;
xx = cumsum(xx);
xx = [0 xx];
for i=1:Nsplits
    JvalRedist{i,2} = Jneg( xx(i)+1:xx(i+1) );
end

J.val   = cell(Nsplits,1);
J.train = cell(Nsplits,1);
for i=1:Nsplits
    J.val{i}   = union(JvalRedist{i,1},JvalRedist{i,2});
    J.train{i} = setdiff(Jall,J.val{i});
end

end