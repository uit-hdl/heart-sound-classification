I = isnan(HSdata.avmeanpg);
[HSdata.murGrade1(I),HSdata.murGrade2(I),...
 HSdata.murGrade3(I),HSdata.murGrade4(I)]

% how many with high murmur grade did not have AS?

% I = HSdata.murGrade1>=2;

%% correlation between AVMEANPG and murmur grade in each position:
clear x U R
I = ~isnan(HSdata.avmeanpg);
[x{1},~,U{1},R{1}] = corrcoef(HSdata.murGrade1(I), HSdata.avmeanpg(I))
[x{2},~,U{2},R{2}] = corrcoef(HSdata.murGrade2(I), HSdata.avmeanpg(I))
[x{3},~,U{3},R{3}] = corrcoef(HSdata.murGrade3(I), HSdata.avmeanpg(I))
[x{4},~,U{4},R{4}] = corrcoef(HSdata.murGrade4(I), HSdata.avmeanpg(I))
[x{5},~,U{5},R{5}] = corrcoef(HSdata.murGradeSum(I), HSdata.avmeanpg(I))
[x{6},~,U{6},R{6}] = corrcoef(HSdata.murGradeMax(I), HSdata.avmeanpg(I))

C = zeros(6,3);
for i=1:6
    C(i,:) = round([U{i}(1,2),x{i}(1,2),R{i}(1,2)]*100,1);
end

array2table(C,'v',["lower" "correlation est." "upper"],...
                'r',["pos. 1" "pos. 2" "pos. 3" "pos. 4" "sum" "max"])



