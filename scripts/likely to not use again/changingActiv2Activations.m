load CVresults_netMurG2AllPos_valStop_overTrain

%% change name of variable
% copy and rename:
CVresults.train.activations = CVresults.train.activ;
CVresults.val.activations   = CVresults.val.activ;

% remove old name:
CVresults.train = rmfield(CVresults.train,'activ');
CVresults.val = rmfield(CVresults.val,'activ');
%% overwrite
save CVresults_netMurG2AllPos_valStop_overTrain CVresults


