
% generate validation and training sets:
load nodes

% remove noisy samples
data0    = HSdataTrainAndVal;
aa       = 4;
gradeThr = 1;
Jnonoise = data0.(sprintf('MURMUR_%gNOISE_REF_T72',aa))==0;
data0           = data0(         Jnonoise,:);
nodes0.loc      = nodes.loc(     Jnonoise,:);
nodes0.state    = nodes.state(   Jnonoise,:);
nodes0.seglines = nodes.seglines(Jnonoise,:);

% DEFINE THE LABELS:
% posCaseStr  = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
% trainOptsScript = 'trainOptsRegress0';
posCaseStr  = 'MRGRADE_T72';
% select training options (Regress for regression, Classif for classification)
trainOptsScript = 'trainOptsClassif0';
Y0      = data0.(posCaseStr)>=1;
N0      = height(data0);
Nsplits = 8;
% number of traning set splits:
Nval = floor(N0/Nsplits);
Jall = 1:N0;

for i=2:Nsplits
% get index for training and validation data
Jval   = (i-1)*Nval+1:i*Nval;
Jtrain = setdiff(Jall,Jval);
Ncomp  = [30,300];                    
trainTime = 60*60;
thrPos    = [];
balanceTrain   = true;
[net,Xtrain,Ytrain] = trainLSTM(data0,Y0,Jtrain,nodes0,aa,...
                        trainOptsScript,trainTime,thrPos,N,balanceTrain,Ncomp);
networks{i} = net;

balanceVal  = false;
[Xval,Yval] = genTrainOrValSet(data0,Y0,Jval,N,aa,balanceVal,...
                            nodes0,Ncomp);
                        
Nloc = 1;
figure
[AUC] = performanceSummaryNeurNet(Xval,Yval,net,...
                                    numel(Jval),N.segPerPCG,true,Nloc);
title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
end

%% compare results with murmur prediction algorithm:
load networksCVnoNoiseDrOutReg
Ypred = predict(networks{1},Xval,'MiniBatchSize',miniBatchSize);
figure
[AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval,Yval=='1',Ypred,...
                                    numel(Jval),N.segPerPCG,true,Nloc);
title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))

                                     