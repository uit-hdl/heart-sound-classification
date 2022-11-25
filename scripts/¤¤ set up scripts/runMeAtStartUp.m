functionFolder = genpath('C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\matlab files\¤¤ functions');  
saveVarFolder = genpath('C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\matlab files\¤¤ saved variables');
setUpScriptsFolder = genpath('C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\matlab files\¤¤ set up scripts');
deepLearningFolder = genpath('C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\matlab files\¤¤ deep learning scripts');
dataFolder = 'C:\Users\perwa\OneDrive - UiT Office 365\Data - Lars Ailo Bongos files';
physionetFolder = genpath('C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\matlab files\¤¤ segmentation');

addpath(functionFolder,saveVarFolder,setUpScriptsFolder,...
        deepLearningFolder,dataFolder,physionetFolder)

clearvars functionFolder saveVarFolder setUpScriptsFolder deepLearningFolder dataFolder physionetFolder
importingAndPreProcessingTUdata
renameAndDefineVariables_TUdata
GenDescriptiveStatisticsTrainingData
GenDescriptiveStatistics
loadHMMpar

% clear variables that are probably not going to be used again:
% initialVars = {'echoData' 'echoDataInfo' 'HMMpar' 'HSdata0' 'HSdata' ...
%     'HSdataTrain' 'Jtest0' 'Jtrain0' 'Jval0' 'murDataInfo' ...
%     'murData' 'stats' 'statsAll' 'initialVars'};

% clear variables that are probably not going to be used again:
initialVars = {'HMMpar' 'HSdata' 'HSdataTrain' 'murDataInfo' 'echoDataInfo' ...
               'Jtest0' 'Jtrain0' 'Jval0' 'stats' 'statsAll' 'initialVars'};
clearvars('-except',initialVars{:})
    