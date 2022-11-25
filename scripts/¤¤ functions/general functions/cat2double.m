function x = cat2double(x)
% convert categorical to double
x = str2double(cellstr(x));
end