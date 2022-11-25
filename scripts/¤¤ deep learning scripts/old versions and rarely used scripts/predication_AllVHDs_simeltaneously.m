% load networksCVnoNoiseDrOutMurG2AllPosValStopOvertrain.mat

% run estSnSpBalancedWay_valStop to get the variable allVHDpredMatrix that
% contains predictions for all VHD's

symptomMatrix = zeros(height(HSdata),4);
symptomMatrix(:,1) = sum([HSdata.ANGINA_T7,HSdata.DYSPNEA_CALMLY_FLAT_T7,...
                        HSdata.CHEST_PAIN_NORMAL_T7],2,'o')>0;
                    

Jall = cell2mat(CVresults.valTot.J);
X = sum(allVHDpredMatrix,2)>0;
X = X(Jall);
Y = [HSdata.ARgrade(Jall)>=4,HSdata.MRgrade(Jall)>=4,...
     HSdata.ASgrade(Jall)>0,HSdata.MSgrade(Jall)>0];

% How many AR are predicted?
condProb(X,Y(:,1))
condProb(~X,~Y(:,1))

% How many MR are predicted?
condProb(X,Y(:,2))
condProb(~X,~Y(:,2))

% How many AS are predicted?
condProb(X,Y(:,3))
condProb(~X,~Y(:,3))

% How many MS are predicted?
condProb(X,Y(:,4))
condProb(~X,~Y(:,4))

% total scores:
condProb(X,sum(Y,2)>0)
condProb(~X,sum(Y,2)==0)
