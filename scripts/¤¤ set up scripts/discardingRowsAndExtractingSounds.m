% in this script I investigate which files have corrupt audio. a total of 5
% observations had corrupt audio, and were removed from the data set
%%  Convert wav. files into time-series
X = cell(height(HSdata0),4);
for i=1:height(HSdata0)
    id = HSdata0.UNIKT_LOPENR(i);
    for j=1:4
        if i==424 || i==588 || i==610 || i==946 || i==1933
            X{i,j} = nan;
        else
            X{i,j} = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,j));
        end
    end
    i
end
%%
JnotFullAudio = [424,588,610,946,1933];
idMissing = HSdata0.UNIKT_LOPENR(JnotFullAudio)
JfullAudio = setdiff(1:height(HSdata),discardedRows)
id = HSdata.UNIKT_LOPENR(discardedRows(2));
HSdata.ARGRADE_T72(HSdata.UNIKT_LOPENR==id)
HSdata.MRGRADE_T72(HSdata.UNIKT_LOPENR==id)
HSdata.ASGRADE_T72(HSdata.UNIKT_LOPENR==id)
HSdata.MSGRADE_T72(HSdata.UNIKT_LOPENR==id)

%%  Convert wav. files into time-series

Xtrain = cell(height(HSdataTrain),4);
for i=1:height(HSdataTrain)
    id = HSdataTrain.UNIKT_LOPENR(i);
    for j=1:4
        if i==351 || i==1022 || i==1120 || i==1452 || i==1692
            Xtrain{i,j} = nan;
        else
            Xtrain{i,j} = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,j));
        end
    end
    i
end


%% investigate the id's that generated errors
% i==351 || i==1022 || i==1120 || i==1452 || i==1692
i = 1692;
id = HSdataTrain.UNIKT_LOPENR(i);
% for j=1:4
%     if 1==2%i==351 || i==1022 || i==1120 || i==1452 || i==1692
%         Xtrain{i,j} = nan;
%     else
%         Xtrain{i,j} = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,j));
%     end
% end

HSdataTrain.ARGRADE_T72(HSdataTrain.UNIKT_LOPENR==id)
HSdataTrain.MRGRADE_T72(HSdataTrain.UNIKT_LOPENR==id)
HSdataTrain.ASGRADE_T72(HSdataTrain.UNIKT_LOPENR==id)
HSdataTrain.MSGRADE_T72(HSdataTrain.UNIKT_LOPENR==id)

HSdataTrain = HSdataTrain(setdiff(1:height(HSdataTrain),[351,1022,1120,1452,1692]),:);