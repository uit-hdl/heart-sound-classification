function TargetArray = getTruthCellArray(dataFrame,IdxArray,TargetVariable,classThr)
% assumes that you have a cell-array of indeces for training or validation
% sets from cross-validation, and want to create a cell array of true
% labels corresponding to those indeces for the target variable called
% targetVariable >= classThr.

%% example
% dataFrame = HSdata;
% IdxArray = CVresults.val.J;
% variable = "murmur";
%%
sz = size(IdxArray);
TargetArray = cell(sz);

for i=1:sz(1)
    for aa=1:sz(2)
        if TargetVariable=="murmur"
            varStr = sprintf('murGrade%g',aa);
        else
            varStr = TargetVariable;
        end
        if isempty(classThr)
            TargetArray{i,aa} = dataFrame.(varStr)(IdxArray{i,aa});
        else
            TargetArray{i,aa} = dataFrame.(varStr)(IdxArray{i,aa})>=classThr;
        end
    end
end

end

