function pval = pValue(x,type)
% gives p-value corresponding to the hypothesis that the elements of x
% comes from a distribution with mean value significantly larger than 0.
if nargin==1
    type = "oneSided";
end

if isrow(x)
    x = x';
end
t = mean(x);
n = height(x);
sigma = std(x);

if type=="oneSided"
    pval = 1 - tcdf(sqrt(n)*t./sigma, n-1);
else
    pval = 2*min([tcdf(sqrt(n)*t./sigma, n-1),1 - tcdf(sqrt(n)*t./sigma, n-1)]);
end

end