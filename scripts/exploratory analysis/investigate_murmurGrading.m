
for i_grade=1:5
    for aa=1:4
        str = sprintf('murGrade%g',aa);
        J{i_grade,aa} = find(HSdata.(str)==i_grade);
        ID{i_grade,aa} = HSdata.id(J{i_grade,aa});
    end
end

%%
aa = 1;
grade = 1;
k = 6;
id = ID{grade,aa}(k);
playHS(id,aa)
close
quickScaleogramPlot(id,aa)