function [phat,phatCi] = calc_SNorSPorProb_CI(Ypred,Ytrue,repr)
% calculates sensitivity or specificity with confidence intervals
% calculated using the binomial distribution.
%% preliminary
if nargin==1
    Ytrue = zeros(size(Ypred))==0;
    repr = [];
elseif nargin==2
    repr = [];
end

if isempty(Ytrue)
    Ytrue = zeros(size(Ypred))==0;
end
%% function code

[phat,phatCi] = getBinomialCI(Ypred==Ytrue);

% if want to display in rounded percentage:
if ~isempty(repr)
    phat = 100*round(phat,repr);
    phatCi = 100*round(phatCi,repr);
end
end