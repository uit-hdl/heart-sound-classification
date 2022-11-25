function [perf,AUC] = balancedPerfEstUsingCVoutput(X,Y,K,balance)
% The elements of X a(i) are validation set ACTIVATIONS contained in a cell
% array generated during K-fold cross-validation. Y is a LOGICAL ground
% TRUTH index-vector. If X is logical, then function returns estimated
% sensitivity, specificity and accuracy. If A is numerial, then function
% outputs estimate of AUC. Rows are SHIFTED as nessecary in order to ensure
% sufficient numbers of POSITIVE cases in each of the K sets after
% reshuffling.
%% example
% X = activ;
% Y = Ytarget;
% K = 8;
%% preliminary
if nargin==2
    % if no K is provided, assume that a metrics are to be computed for the
    % combined validation predictions.
    K = 1;
end
%% code
if K==1
    J.val{1} = 1:numel(Y);
else
    J = distributePosCasesEvenly(Y,K);
end

XcellBal = cell(K,1);
YcellBal = cell(K,1);
% reshuffle rows to ensure sufficient numbers of positive cases in all sets:
for i=1:K
    XcellBal{i} = X(J.val{i});
    YcellBal{i} = Y(J.val{i});
end

if isnumeric(X)
    % estimate AUC:
    AUC = zeros(K,1);
    for i=1:K
        [~,~,~,AUC(i)] = perfcurve(YcellBal{i},XcellBal{i}, true);
    end
    perf = computeCImeanEst(AUC, "2","tdist");
    
elseif islogical(X)
    % estimate SN, SP and AC:
    perf = cell(K,1);
    NcorrPred = zeros(3,1);
    Ntrials   = zeros(3,1);
    for i=1:K
        predictions = XcellBal{i};
        Ytrue       = YcellBal{i}; 
        % SN:
        NcorrPred(1) = sum(and(predictions,Ytrue));
        Ntrials(1)   = sum(Ytrue);
        % SP:
        NcorrPred(2) = sum(and(~predictions,~Ytrue));
        Ntrials(2)   = sum(~Ytrue);
        % AC:
        NcorrPred(3) = sum(predictions==Ytrue);
        Ntrials(3)   = numel(Ytrue);

        [phat,phatCi] = binofit(NcorrPred,Ntrials);
        perf{i} = [phat';phatCi'];
    end
    
end

end




