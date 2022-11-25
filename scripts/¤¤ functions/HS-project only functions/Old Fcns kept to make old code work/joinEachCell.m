function Z = joinEachCell(S)
% takes a structure S that conttains two cell arrays, S.activations and S.J
% that is obtained during cross-validation, and merges them into a single
% cell-array of two-coulmn arrays. Each cell array is assumed to be of same
% size N-by-4, and cell i,j in S.activations is assumed to have the same
% dimensions as cell i,j in S.J
%% example
% S = CVresults.val;
%%
Nsplits = height(S.J);
Z = cell(Nsplits,4);
for i=1:Nsplits
    for aa=1:4
        Z{i,aa} = [S.J{i,aa},S.activations{i,aa}];
    end
end

end