function diseaseIndex = disease2index(diseaseName)
% convenience function that maps name of disease to an index value. Mostly
% used in loops.
if diseaseName=="AR"
    diseaseIndex = 1;
elseif diseaseName=="MR"
    diseaseIndex = 2;
elseif diseaseName=="AS"
    diseaseIndex = 3;
elseif diseaseName=="MS"
    diseaseIndex = 4;
elseif diseaseName=="avmpg"
    diseaseIndex = 1;
elseif diseaseName=="ASpg"
    diseaseIndex = 1;
elseif diseaseName=="ASnew"
    diseaseIndex = 1;
elseif contains(diseaseName,"VHD")
    diseaseIndex = 5;
elseif diseaseName=="costumTarget"
    diseaseIndex = 6;
end

end