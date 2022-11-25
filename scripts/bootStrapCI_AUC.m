


myFunc = @(x1,x2) getAUC(x1,x2);
% myFunc = @(x1,x2) condProb(x2>=.3,x1);
nIterations = 100;
% compute the CI of correlation between murmur grade and
% aortic-valve-pressure-gradient:
CI = bootci(nIterations,...
    {myFunc, Ytarget, activ})