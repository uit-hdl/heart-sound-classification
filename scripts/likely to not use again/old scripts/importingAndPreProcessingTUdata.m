% Convert variable names 
opts = detectImportOptions('data_lars_27Sep2021.csv');
opts = detectImportOptions('variables-sounds-original.xlsx');
opts.DataLines = [1,1];

varNamesCell = readtable('data_lars_27Sep2021.csv',opts);
varNamesCell = readtable('variables-sounds-original.xlsx',opts);
varNamesCell = table2cell(varNamesCell);
opts.DataLines = [2,inf];
opts.VariableNames = varNamesCell;

opts = setvartype(opts,varNamesCell, 'string');
TUdata = readtable('data_lars_27Sep2021.csv',opts);

varNamesStr = string(varNamesCell);
%% Extract the rows with murmur-data
Imurmur = zeros(1,size(TUdata,1));
for i=1:size(TUdata,1)
    Imurmur(i) = length( str2num(TUdata.MURMUR_1NORMAD_T72(i)) );
end
Jmurmur = find(Imurmur)';
HSdata0 = TUdata(Jmurmur,:); % original data set with HS Ausc rows

%% Remove Nonsensical Data-rows
% some of the rows do not make any sense! collect the ones that do make
% sense:
cleanId = zeros(1,size(HSdata0,1));
for i=1:size(HSdata0,1)
    cleanId(i) = length( str2num(HSdata0.idrand(i) ));  %#ok<*ST2NM>
end

HSdata = HSdata0(logical(cleanId),:);
%% Extract and describe murmur data and save as separate table
murData = HSdata;
genMurmDataInfoStructure
HSdata = murData;
murData = murData(:,murDataInfo.allMurData.col);
%% Extract and describe ECHO-data and save as separate table
echoData = HSdata;
genEchoDataInfoStructure
echoData = HSdata(:,echoDataInfo.allEchoData.col);

%% shift column indeces to match the Murmur-data-table
fn = fieldnames(murDataInfo);
Ntemp = width(HSdata)-width(HSdata0);
for k=2:numel(fn)
    if k<22 || k>24
        if isfield(murDataInfo.(fn{k}),'newData')
            murTableInd.(fn{k}) = murDataInfo.(fn{k}).col - (257-146) - 40;
        else
            murTableInd.(fn{k}) = murDataInfo.(fn{k}).col - 40;
        end
    end
end
%% shift column indeces to match the Echo-data-table
fn = fieldnames(echoDataInfo);
for k=2:numel(fn)
    echoTableInd.(fn{k}) = echoDataInfo.(fn{k}).col - 146;
end
%% Find for which ID's echo was performed
echoRows = not(ismissing(echoData.ECHO_CHKIN_DATE_T72));
echoRows = logical(echoRows);
echoDataFull = echoData;
echoDataRed  = echoData(echoRows,:);
murDatafull = murData;
murDataRed  = murData(echoRows,:);
HSdataRed  = HSdata(echoRows,:);
%% convert explanatory variables from string to double
HSdataRed = convertvars(HSdataRed,varNamesCell(1:40),'double');
%% convert continuous echo raw data into doubles
HSdataRed = convertvars(HSdataRed,varNamesCell([166,194:257]),'double');
%% For VHD grades, replace nan with 0
cols = echoTableInd.stenosisAndRegurgGrade;
Inan = isnan(echoDataRed{:,cols});

echoDataRed{Inan(:,1),cols(1)} = 0;
echoDataRed{Inan(:,2),cols(2)} = 0;
echoDataRed{Inan(:,3),cols(3)} = 0;
echoDataRed{Inan(:,4),cols(4)} = 0;

cols = echoDataInfo.stenosisAndRegurgGrade.col;
Inan = isnan(HSdataRed{:,cols});

HSdataRed{Inan(:,1),cols(1)} = 0;
HSdataRed{Inan(:,2),cols(2)} = 0;
HSdataRed{Inan(:,3),cols(3)} = 0;
HSdataRed{Inan(:,4),cols(4)} = 0;

cols = 194:256;
Inan = isnan(HSdataRed{:,cols});
%% split into validation set and training set
load('Jtrain.mat')
load('Jval.mat')
HSdataTrain = HSdataRed(Jtrain,:);
HSdataVal   = HSdataRed(Jval,:);





