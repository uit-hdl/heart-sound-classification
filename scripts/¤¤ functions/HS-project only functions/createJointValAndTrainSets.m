function new = createJointValAndTrainSets(CVresults)
% trainSets and Valsets are Nsplits-by-4 cells that contains linear indeces
% for training and validation sets obtained during cross validation for
% each  of the 4 auscultation positions. The indeces of the CV sets may not
% agree exactly since some recording, and to fix this, indeces from the training sets are added
% to individual validation sets and subtracted from the corresponding
% training sets in order to create joint sets.
%% define more convenient notation
JtrainSets = CVresults.train.J;
JvalSets   = CVresults.val.J;
activTrain = CVresults.train.activations;
activVal   = CVresults.val.activations;
%%
Nsplits = height(CVresults.train.J);
old.trainMat = cell(Nsplits,4);
old.valMat   = cell(Nsplits,4);

for i=1:Nsplits
    for aa=1:4
        old.trainMat{i,aa} = [JtrainSets{i,aa},activTrain{i,aa}];
        old.valMat{i,aa}   = [JvalSets{i,aa}  ,activVal{i,aa}];
    end
end
new = old;

for i=1:Nsplits
    
    JunionVal   = unionIterated(JvalSets(i,:));
    JunionTrain = unionIterated(JtrainSets(i,:));
    % define joint validation and training sets:
    coreValSet = intersectionIterated(JvalSets(i,:));
    coreValIndex = median(coreValSet);
    
    % the mixed set contains all indeces with both training and validation
    % predictions.
    mixedSet   = intersect(JunionVal,JunionTrain);
    mixedSetL  = mixedSet(mixedSet < coreValIndex);
    mixedSetR  = mixedSet(mixedSet > coreValIndex);
    
    for aa=1:4
        % index of what gets added to validation and removed from training:
        JvalAdd = intersect(JtrainSets{i,aa},mixedSetR);
        % index of what gets removed from validation and added to training:
        JvalSub = intersect(JvalSets{i,aa},mixedSetL);
        
        % *** add to validation and subtract from training set ***
        Iadd = findInd(JtrainSets{i,aa},JvalAdd);
        shiftRowsR = old.trainMat{i,aa}(Iadd,:);
        % add rows from right to validation set:
        new.valMat{i,aa}   = [old.valMat{i,aa};shiftRowsR];
        % remove rows from training set:
        new.trainMat{i,aa} = old.trainMat{i,aa}(~Iadd,:);
        
        % *** add to training and subtract from validation set ***
        Isub = findInd(new.valMat{i,aa}(:,1),JvalSub);
        shiftRowsL = new.valMat{i,aa}(Isub,:);
        % add rows from left to training set:
        new.trainMat{i,aa} = [shiftRowsL; new.trainMat{i,aa}];
        % remove rows from validation set:
        new.valMat{i,aa} = new.valMat{i,aa}(~Isub,:);
    end
end

end















%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD VERSION 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% function mod = createJointValAndTrainSets(JtrainSets,JvalSets,activTrain,activVal)
%%
% JtrainSets = CVresults.Jtrain;
% JvalSets = S.indiv.J;
% activTrain = CVresults.activations;
% activVal = S.indiv.activation;
%%
% Nsplits = height(JtrainSets);
% mod.JtrainSets = cell(Nsplits,4);
% mod.JvalSets   = cell(Nsplits,4);
% mod.activTrain = cell(Nsplits,4);
% mod.activVal   = cell(Nsplits,4);
% 
% JvalSpans = cell(Nsplits,1);
% for i=1:Nsplits
%     
%     JunionValSets = unionIterated(JvalSets(i,:));
%     JvalSetsSpan  = min(JunionValSets):max(JunionValSets);
%     JvalSpans{i} = JvalSetsSpan;
%     
%     % *** add to validation set from above ***
%     for aa=1:4
%         % get shift indices:
%         shiftInd  = intersect(JtrainSets{i,aa},JvalSetsSpan);
%         if i>1
%             shiftInd = setdiff(shiftInd,JvalSpans{i-1});
%         end
%         % get position of shift indeces within the training set:
%         Ishift = findInd(JtrainSets{i,aa},shiftInd);
%         % get shift activations:
%         shiftActiv  = activTrain{i,aa}(Ishift);
%         % add activations to validation set:
%         mod.activVal{i,aa} = [activVal{i,aa};shiftActiv];
%         % subtract activations from validation set:
%         mod.activTrain{i,aa} = activTrain{i,aa}(~Ishift);
%         % add shift indices to validation set:
%         mod.JvalSets{i,aa}   = union(JvalSets{i,aa},shiftInd);
%         % subtract shift indeces from training set:
%         mod.JtrainSets{i,aa} = setdiff(JtrainSets{i,aa},shiftInd);
%     end
%     
%     % *** subtract from validation set from below ***
%     for aa=1:4
%         
%     end
%         
% end

%%%%%%%%%%%%%% OLD VERSION 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% %%
% Nsplits = height(JtrainSets);
% mod.JtrainSets = cell(Nsplits,4);
% mod.JvalSets   = cell(Nsplits,4);
% mod.activTrain = cell(Nsplits,4);
% mod.activVal   = cell(Nsplits,4);
% 
% for i=1:Nsplits
% 
%     JunionValSets = unionIterated(JvalSets(i,:));
%     JvalSetsSpan  = min(JunionValSets):max(JunionValSets); 
%     for aa=1:4
%         % get shift indices:
%         shiftInd  = intersect(JtrainSets{i,aa},JvalSetsSpan);
%         % separate shift indeces into those that gets added from below and
%         % above respectively:
%         shiftIndBelow = shiftInd(shiftInd<min(JvalSets{i,aa}));
%         shiftIndAbove = shiftInd(shiftInd>max(JvalSets{i,aa}));
%         % get position of shift indeces within the training set:
%         IshiftBelow = findInd(JtrainSets{i,aa},shiftIndBelow);
%         IshiftAbove = findInd(JtrainSets{i,aa},shiftIndAbove);
%         % get shift activations:
%         shiftActivBelow  = activTrain{i,aa}(IshiftBelow);
%         shiftActivAbove  = activTrain{i,aa}(IshiftAbove);
%         
%         mod.activVal{i,aa} = [shiftActivBelow;activVal{i,aa};shiftActivAbove];
%         mod.activTrain{i,aa} = activTrain{i,aa}(~or(IshiftBelow,IshiftAbove));
%         % subtract shift indeces from training set:
%         mod.JtrainSets{i,aa} = setdiff(JtrainSets{i,aa},shiftInd);
%         % add shift indices to validation set:
%         mod.JvalSets{i,aa}   = union(JvalSets{i,aa},shiftInd);
%     end
% end
