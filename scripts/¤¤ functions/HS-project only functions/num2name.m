function fieldName = num2name(k,preString)
% convert an integer to a string that can be used as an acceptable field
% name of the form "g1", as in grade 1. Makes life ever so slightly easier.

% ** EXAMPLE OF USAGE **
% Say you have a field F, and you want it to have a field called g2. then
% run num2name(2), and the function returns 'g2'. Convenient to use in
% loops.

if nargin==1
    preString = 'g'; % g for grade
end
fieldName = strcat(preString,num2str(k));
end

