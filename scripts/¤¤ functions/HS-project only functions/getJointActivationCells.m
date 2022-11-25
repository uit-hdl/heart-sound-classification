function jointActivCell = getJointActivationCells(indivActivCell)
% indivActivCell contains two 8-by-4 cell arrays, one for the validation sets
% and one for the training sets, and is assumed to be the output of the
% function createJointValAndTrain, so that training predictions have been
% shifted to validation sets for prediction-padding. Outputs two 1 by 8
% cell arrays, one for validation-sets and one for training-sets, each of
% which contains indeces and max-activations for the corresponding set.

%% example
% indivActivCell = new;
% get range
%%
Nsplits = height(indivActivCell.trainMat);
jointActivCell.train = cell(Nsplits,1);
jointActivCell.val   = cell(Nsplits,1);
Nsplits = height(indivActivCell.trainMat);

for k=1:2
    if k==1
        CVcellJoint = indivActivCell.valMat;
    else
        CVcellJoint = indivActivCell.trainMat;
    end

    for i=1:Nsplits
    cellRow_i = CVcellJoint(i,:);
    % get indeces for observations with at least one prediction:
    Jjoint = cellRow_i{1}(:,1);
    for aa=2:4
        Jjoint = union(Jjoint, cellRow_i{aa}(:,1));
    end

    % get span of observations:
    S = zeros(4,2);
    for aa=1:4
        S(aa,:) = getSpan(cellRow_i{aa}(:,1));
    end
    S = [min(S(:,1)),max(S(:,2))];
    activMatrix = zeros(range(S)+1,4);

    for aa=1:4
        activMatrix(cellRow_i{aa}(:,1)-S(1)+1,aa) = cellRow_i{aa}(:,2);
    end

    % remove rows that contain only noise recordings (all 0s):
    activMatrix = activMatrix(sum(activMatrix,2)>0,:);
    % let maximum activation represent the joint activation:
    jointActivVec = max(activMatrix,[],2);
    
    if k==1
        jointActivCell.val{i}   = [Jjoint,jointActivVec];
    else
        jointActivCell.train{i} = [Jjoint,jointActivVec];
    end
    
    end


end

end
