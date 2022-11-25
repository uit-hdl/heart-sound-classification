load CVresults_netMurRegAllPos_valStop_overTrain

% *** how many murmurs detected where innocent? ***
allActivs = reshape(CVresults.val.activations,[8*4,1])
allJs     = reshape(CVresults.val.J,[8*4,1])
allActivs = cell2mat(allActivs);
allJs = cell2mat(allJs);

allYtargets = (HSdata.ASgrade(allJs)>0) + (HSdata.MSgrade(allJs)>0) + ...
              (HSdata.ARgrade(allJs)>=2) + (HSdata.MRgrade(allJs)>=2);
allYtargets = (allYtargets)>0

% what proportion was innocent murmur?
condProb(~allYtargets,allActivs>=2)
