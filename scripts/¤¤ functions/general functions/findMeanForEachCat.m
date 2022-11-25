function [meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(x,y)
% ** description **
% x is an integer valued vector representing categories, and y is numerical
% vector with values corresponding to the categories in x. Returns mean
% values for each category, ordered from smallest to largest element in x.
% Also returns a cell array with index vectors for each integer category;
% ;a cell array with the y-values for each category; and a cell array
% containing confidence intervals for each category. Takes mean of the
% Ignores nan-values.

% ** used for **
% example input:
%  x = [1   1  2  3  2  2 1   3];
%  y = [.1 .2 .5 .5 .1 .1 1.0 1.5]
% ...
%%
% find integers in x
cat = unique(x);
nCat = length(cat);
Icat = cell(1,nCat);
yCat = cell(1,nCat);
CI   = zeros(2,nCat);
meanVals = zeros(1,nCat);

for i=1:nCat
    Icat{i} = (x==cat(i));
    yCat{i} = y(Icat{i});
    meanVals(i) = mean(yCat{i},'omitnan');
    CI(:,i) = computeCImeanEst(yCat{i});
end


end