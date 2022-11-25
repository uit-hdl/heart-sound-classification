%%
g =@(x,a,b,c) x.^(a*exp(-b*x));
a = 1;
b = 0.000;
clf
[meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(DATA.maxMeanMurGrade,DATA.ARGRADE_T72);
plot(DATA.maxMeanMurGrade + randn(height(DATA),1)*0.1, g(DATA.ARGRADE_T72,a,b),'.')
hold on
scatter(cat,g(meanVals,a,b),[],'k','filled')
%%

lm1 = fitlm(DATA,'AR ~ maxMeanMurGrade + AGE_T7 + SEX_T7');

lm2 = fitlm(X,Y);