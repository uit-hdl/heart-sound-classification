function y = myunique(x)
% costume made version of "unique", which only returns one nan-value.
y = unique(x);
if iscategorical(y)
    if any(ismissing(y))
        y(ismissing(y)) = []; % remove all nans
        y(end+1) = missing; % add the unique one.
    end
else
    if any(isnan(y))
        y(isnan(y)) = []; % remove all nans
        y(end+1) = NaN; % add the unique one.
    end
end
end