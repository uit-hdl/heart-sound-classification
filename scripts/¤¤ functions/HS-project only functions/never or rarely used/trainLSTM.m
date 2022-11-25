function [net,Xtrain,Ytrain] = trainLSTM(dataFrame,Y0,Jtrain,nodes,aa,...
                                trainOpts,trainTime,thrPos,N,balance,Ncomp)
%% to use default settings only, use:
% trainLSTM(dataTrain,Ytrain,aa,thrPos,nodes,Jtrain)

% inputs:
% dataFrame = data frame from which the training set is obtained.
% aa     = auscultation location to train on
% Jtrain = index vector
% thrPos = threshold that seperates positive and negative cases
% nodes  = field that contains segmentation information. Assumed to
%          correspond to the data in dataFrame.
% Y0     = output data/labels. can be numeric of logical
% N      = field that contains information about how much to downsample
%           timeseries and segment extraction.
% balance = specifies whether or not to resample to balance the dataset and
%           ensure equal number of positive and negative cases.
% learnOpt = options for learning.
%% default values
if isempty(N)
    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;
end
if isempty(Jtrain)
    Jtrain = 1:height(dataTrain);
end
if isempty(thrPos)
    thrPos = 1;
end
if isempty(balance)
    balance = true;
end
if isempty(Ncomp)
    Ncomp = [30,200];
end

%% Get training set:
[Xtrain,Ytrain] = genTrainOrValSet(dataFrame,Y0,Jtrain,N,aa,...
                                  balance,nodes,Ncomp,thrPos);
                        
%% run training option script to define architecture and get training settings
run(trainOpts);
%% Train:
net = trainNetwork(Xtrain,Ytrain,layers,options);

end