function u_opt = getOptimalThr(X,Y,T,minSn,minSp)
% takes sensitivity and specificity for different thresholds stored in X
% and Y (the output from performanceSummaryNeurNet) and finds the threshold
% u_opt that maximizes the sum of the sensitivity and specificity, subject
% to the condition that the sensitivity must be <= minSensitivity.
if nargin==3
    minSn = 0;
    minSp = 0;
elseif nargin==4
    minSp = 0;
end

if isempty(X) % make first argument empty if you want to compute X,Y and T
    Ytarget = Y;
    activ = T;
    [~,X,Y,T,~] = performanceSummaryNeurNet([],Ytarget,activ,...
                                            [],[],false);
    spec = 1-X.whole;
    sens = Y.whole;
    T = T.whole;
elseif isstruct(X)
    spec = 1-X.whole;
    sens = Y.whole;
    T = T.whole;
else
    spec = 1 - X;
    sens = Y;
end

% find the activation-thresholds that satisfies minimum requirements:
Ifeasible = and(sens>=minSn, spec>=minSp);

if sum(Ifeasible)==0
    warning('no threshold satisfies criteria')
    Ifeasible = sens>=minSn;
end

[~,i_optimal] = max((sens + spec).*Ifeasible);
% get decision optimal threshold:
u_opt = T(i_optimal);

end