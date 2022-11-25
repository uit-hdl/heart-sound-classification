function fileName = getNameOfWav(set,k)
% Convenience function for getting the names of the wav. files in the
% physionet 2016 training set.Possible sets are a,b,c,d,e and f.

% ¤¤¤ sizes of datasets ¤¤¤ 
% A 65
% b 
if set=='e'
    s = '0';
else
    s = '';
end

if k<10
    fileName = sprintf('%s%s000%g.wav',set,s,k);
elseif k>=10 && k<100
    fileName = sprintf('%s%s00%g.wav',set,s,k);
else
    fileName = sprintf('%s%s0%g.wav',set,k);
end

end