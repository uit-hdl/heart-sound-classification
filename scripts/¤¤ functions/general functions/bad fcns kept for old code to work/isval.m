function y = isval(x)
% finds non-nan values
if isnumeric(x)
    y = ~isnan(x);
else
    y = ~ismissing(x);
end
end