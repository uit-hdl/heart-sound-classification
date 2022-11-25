% This is the script were the data is divided into training, validation and
% test set. The model developement data is train + val = 80% of the data.
% Variables are stored in dataSplitID.mat and dataSplitJ.mat

% 90% if data is set to model developement data (train+val), 10% is test
% data:
Ntot   = height(HSdata);
Ntrain = floor(height(HSdata)*0.80);
Nval = floor(height(HSdata)*0.10);
Ntest = Ntot - Ntrain - Nval;

% get linear row indeces
Jfull  = 1:Ntot;
Jtrain = randsample(Ntot,Ntrain);
Jval   = randsample(setdiff(Jfull,Jtrain),Nval);
Jtest   = setdiff(Jfull,union(Jtrain,Jval));

% get participant id's:
IDtrain = HSdata.UNIKT_LOPENR(Jtrain);
IDval   = HSdata.UNIKT_LOPENR(Jval);
IDtest  = HSdata.UNIKT_LOPENR(Jtest);
%%
HSdataTrain = HSdata(Jtrain,:);
HSdataVal   = HSdata(Jval,:);
HSdataTest  = HSdata(Jtrain,:);
%% save data
save('dataSplitJ.mat','Jtrain','Jval','Jtest')
save("dataSplitID.mat",'IDtrain','IDval','IDtest')



