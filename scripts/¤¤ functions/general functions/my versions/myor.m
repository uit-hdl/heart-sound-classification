function Ior = myor(indexArray,arg2,nanTreatment)
% my version of the function "or", which accepts as argument a cell array
% of logical index vectos. 

% nanTreatment is an optional argument that determines how to treat entries
% where at least one of the vectors has a nan value. The default is for
% both values to set as nan only when both are nan ("or").

% *** Possible options ***
% "and" --> put nan only when both have nan.
% "or"  --> put nan if nan in one or more of the vectors

if nargin<3
    nanTreatment = "or";
end

if iscell(indexArray)
    Ior = indexArray{1}==1;
    if length(indexArray)>1
        for i=2:length(indexArray)
            Ior = or(Ior==1,indexArray{i}==1);
        end
    end
    
else 
    
    arg1 = indexArray;
    
    Inan = nanVec(arg1,arg2,nanTreatment);    
    Ior = or(indexArray==1,arg2==1) + Inan;

end

end