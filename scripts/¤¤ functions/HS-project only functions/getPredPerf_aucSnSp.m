function [T,M] = getPredPerf_aucSnSp(Ypred,Ytrue,AUCmat,nsigdig,calcAcc)
% computes sn, sp, and AUC (and accuracy if requested) along with
% confidence intervals, and collects these in a table. Inputting nsigdig
% informs the function that the output should be displayed in percentage
% format, and should be rounded to nsigdig significant digits. Ytrue and
% Ypred should be in reduced format, and AUC mat should contain the AUC
% values produced during cross validation. 
if nargin==3
    nsigdig = [];
    calcAcc = false;
elseif nargin==4
    calcAcc = false;
end

[AUCci,AUCciWidth,AUC] = computeCImeanEst(AUCmat);
computeCI = true;
[sn,sn_ci] = condProb(Ypred,Ytrue,     computeCI,[],"binomial");
[sp,sp_ci] = condProb(~Ypred,~Ytrue,   computeCI,[],"binomial");
[ac,ac_ci] = condProb(Ypred==Ytrue ,[],computeCI,[],"binomial");

M = [[sn       sp       ac        AUC];...
     [sn_ci(1) sp_ci(1) ac_ci(1)  AUCci(1)];...
     [sn_ci(2) sp_ci(2) ac_ci(2)  AUCci(2)];...
     [nan      nan      nan       AUCciWidth]];

if not(calcAcc)
    M = M(:,[1,2,4]);
end
 
if ~isempty(nsigdig)
    M = round(100*M,nsigdig);
end

T = array2table(M,'v',{'sn','sp','AUC'},'r',{'est.','ciLower','ciUpper','ciWidth'});

end
