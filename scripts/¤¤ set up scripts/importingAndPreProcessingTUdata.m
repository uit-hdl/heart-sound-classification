% Convert variable names 
opts = detectImportOptions('variables-sounds-original.xlsx');

TUdata = readtable('variables-sounds-original.xlsx',opts);
varNamesCell = opts.SelectedVariableNames;

% convert MSGRADE from string to double:
TUdata.MSGRADE_T72 = str2double(TUdata.MSGRADE_T72);
%% Extract the rows with murmur-data
HSdata = TUdata(~isnan(TUdata.MURMUR_1NORMAD_T72),:); % original data set with HS Ausc rows

%% Correct for the fact that ausc-location ordering was changed during data collection:
% correct for the fact that the ordering of the auscultation areas was
% reversed after week 34:
load('idAuscOrder.mat','idAuscOrder')
HSdata.UNIKT_LOPENR = str2double(HSdata.UNIKT_LOPENR);
correctOrderOfMurmurData
%% Extract and describe murmur data and save as separate table
murData = HSdata;
genMurmDataInfoStructure
HSdata  = murData;
murData = murData(:,murDataInfo.allMurData.col);
%% Extract and describe ECHO-data and save as separate table
HSdata(HSdata.UNIKT_LOPENR==10492521,:) = []; % Ej fullstendig echo information, lettest å bare plukke bort

HSdata.ARGRADE_T72(isnan(HSdata.ARGRADE_T72)) = 0;
HSdata.MRGRADE_T72(isnan(HSdata.MRGRADE_T72)) = 0;
HSdata.ASGRADE_T72(isnan(HSdata.ASGRADE_T72)) = 0;
HSdata.MSGRADE_T72(isnan(HSdata.MSGRADE_T72)) = 0;

HSdata.MSGRADE_T72 = convertCharsToStrings(HSdata.MSGRADE_T72);
%% extract echo-data
echoData = HSdata; %#ok<*NASGU>
genEchoDataInfoStructure
echoData = HSdata(:,[1,echoDataInfo.allEchoData.col]); % the one is to keep the id-column
%% Remove rows that have incomplete audio data
load('usableRows.mat','usableRows')
HSdata0 = HSdata;
HSdata = HSdata0(usableRows,:);
% discardedRows = setdiff(1:height(HSdata0),usableRows);
% Nusable rows = 2124
% Ntotal = 2129
% Ndiscarded = 5

%% split into validation set and training set
% note that cross-validation set = original training set + original val-set
load('dataSplitJ.mat')
load('dataSplitID.mat')
HSdataTrain = HSdata(union(Jtrain0,Jval0),:);

