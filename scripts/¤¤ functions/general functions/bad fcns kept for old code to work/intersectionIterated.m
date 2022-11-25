function Jintersection = intersectionIterated(Jcell)
% convenience function that allows you to take the union of arbitrary
% numbers of index vectors instead of just two.

n = max(size(Jcell));

if isnumeric(Jcell{1})
    Jintersection = Jcell{1};
    for i=2:n
        Jintersection = intersect(Jintersection,Jcell{i});
    end
else
    Jintersection = find(Jcell{1});
    for i=2:n
        Jintersection = intersect(Jintersection,find(Jcell{i}));
    end
end

end