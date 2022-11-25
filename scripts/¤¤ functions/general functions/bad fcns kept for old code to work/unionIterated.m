function Junion = unionIterated(Jcell,outputType)
% convenience function that allows you to take the union of arbitrary
% numbers of index vectors instead of just two.

if nargin==1
    outputType = "linear";
end

n = max(size(Jcell));

if isnumeric(Jcell{1})
    Junion = Jcell{1};
    for i=2:n
        Junion = union(Junion,Jcell{i});
    end
elseif islogical(Jcell{1})
    Junion = Jcell{1};
    for i=2:n
        Junion = Junion + Jcell{i};
    end
    Junion = Junion>0;
    if outputType=="linear"
        Junion = find(Junion);
    end
end

end