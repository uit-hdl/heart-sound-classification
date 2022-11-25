DATA = HSdataTrain(not(isnan(HSdataTrain.AVMAXPG_T72)),:);
DATA = DATA(~isnan(DATA.ARPRESENCE_T72),:);
dataNames = {'maxMeanMurGrade','AGE_T7','SEX_T7','DIABETES_T7'...
    'DYSPNEA_FAST_UPHILL_T7','AVMAXPG_T72'};

 
clear X
% DATA.ASPRESENCE_T72 = categorical(DATA.ASPRESENCE_T72);
X(:,1) = DATA.(dataNames{1});
X(:,2) = DATA.(dataNames{2});
X(:,3) = categorical(DATA.(dataNames{3}));
X(:,4) = categorical(DATA.(dataNames{4}));
X(:,5) = DATA.(dataNames{5});
X(:,6) = DATA.(dataNames{6});
Nvars = size(X,2);
Y = DATA.AVMAXPG_T72.^.2;
    %%
clf
SPS = findSubplotSize(Nvars);
for i=1:Nvars
    subplot(SPS(1),SPS(2),i)
    plot(X(:,i),Y,'o')
    title(sprintf('%s',dataNames{i}))
end
%% MAKE LINEAR MODEL FOR PREDICTING AS PRESENCE BASED ON MURMUR GRADE 
g =@(x,a,b,c) x.^(a*exp(-b*x));
a = 0.25;
b = 0.08;
clf
[meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(DATA.maxMeanMurGrade,DATA.ASGRADE_T72)
plot(DATA.maxMeanMurGrade,g(DATA.ASGRADE_T72,a,b),'.')
hold on
scatter(cat,g(meanVals,a,b),[],'k','filled')
%%
DATA = HSdataTrain(not(isnan(HSdataTrain.AVMAXPG_T72)),:);
DATA = DATA(~isnan(DATA.ARPRESENCE_T72),:);
DATA.DYSPNEA_FAST_UPHILL_T7 = categorical(DATA.DYSPNEA_FAST_UPHILL_T7);
DATA.HIGH_BLOOD_PRESSURE_T7 = categorical(DATA.HIGH_BLOOD_PRESSURE_T7);
DATA.SMOKE_DAILY_Q2_T7 = categorical(DATA.SMOKE_DAILY_Q2_T7);
DATA.DIABETES_T7 = categorical(DATA.DIABETES_T7);
DATA.SEX_T7 = categorical(DATA.SEX_T7);

lm1 = fitlm(DATA,'AVMAXPG_T72 ~ maxMeanMurGrade + AGE_T7 + SEX_T7 + BMI_T7 + DYSPNEA_FAST_UPHILL_T7 + DIABETES_T7 + SMOKE_DAILY_Q2_T7');
DATA.AVMAXPG_T72 = g(DATA.ASGRADE_T72,a,b);
lm2 = fitlm(DATA,'AVMAXPG_T72 ~ maxMeanMurGrade + AGE_T7 + SEX_T7');

lm3 = fitlm(X,Y);
%%
% I = lm2.Fitted>3;
% s1(1) = sum(and(DATA.ARPRESENCE_T72,I))/sum(DATA.ARPRESENCE_T72)
% s1(2) = sum(and(~DATA.ARPRESENCE_T72,~I))/sum(~DATA.ARPRESENCE_T72)
mean(lm2.Fitted)
I = lm2.Fitted>mean(lm2.Fitted)*(1+0.15);
trueLab = DATA.ASPRESENCE_T72;
s2(1) = sum(and(trueLab,I))/sum(trueLab); %#ok<*NOPTS>
s2(2) = sum(and(~trueLab,~I))/sum(~trueLab);
posPredVal = sum(and(DATA.ASPRESENCE_T72,Ipred))/sum(Ipred);



%% MAKE LINEAR MODEL FOR PREDICTING MS PRESENCE BASED ON MURMUR GRADE %%%%%
g =@(x,a,b,c) x.^(a*exp(-b*x));
a = 0.25;
b = 0.08;
clf
[meanVals,Icat,yCat,CI,cat] = findMeanForEachCat(DATA.maxMeanMurGrade,DATA.MSGRADE_T72)
plot(DATA.maxMeanMurGrade,g(DATA.MSGRADE_T72,a,b),'.')
hold on
scatter(cat,g(meanVals,a,b),[],'k','filled')
%%
DATA = HSdataTrain(not(isnan(HSdataTrain.AVMAXPG_T72)),:);
DATA = DATA(~isnan(DATA.ARPRESENCE_T72),:);
DATA.DYSPNEA_FAST_UPHILL_T7 = categorical(DATA.DYSPNEA_FAST_UPHILL_T7);
DATA.HIGH_BLOOD_PRESSURE_T7 = categorical(DATA.HIGH_BLOOD_PRESSURE_T7);
DATA.SMOKE_DAILY_Q2_T7 = categorical(DATA.SMOKE_DAILY_Q2_T7);
DATA.DIABETES_T7 = categorical(DATA.DIABETES_T7);
DATA.SEX_T7 = categorical(DATA.SEX_T7);

lm1 = fitlm(DATA,'AVMAXPG_T72 ~ maxMeanMurGrade + AGE_T7 + SEX_T7 + BMI_T7 + DYSPNEA_FAST_UPHILL_T7 + DIABETES_T7 + SMOKE_DAILY_Q2_T7');
DATA.AVMAXPG_T72 = g(DATA.MSGRADE_T72,a,b);
lm2 = fitlm(DATA,'MSGRADE_T72 ~ maxMeanMurGrade + AGE_T7 + SEX_T7');

lm3 = fitlm(X,Y);
%%
% I = lm2.Fitted>3;
% s1(1) = sum(and(DATA.ARPRESENCE_T72,I))/sum(DATA.ARPRESENCE_T72)
% s1(2) = sum(and(~DATA.ARPRESENCE_T72,~I))/sum(~DATA.ARPRESENCE_T72)
mean(lm2.Fitted)
Ipred = lm2.Fitted > mean(lm2.Fitted)*(1-0.3);
trueLab = DATA.MSPRESENCE_T72;
s2(1) = sum(and(trueLab,Ipred))/sum(trueLab) %#ok<*NOPTS>
s2(2) = sum(and(~trueLab,~Ipred))/sum(~trueLab)
posPredVal = sum(and(DATA.MSPRESENCE_T72,Ipred))/sum(Ipred);






%%
pThr = mean(lm2.Fitted)*0.1:.0001:mean(lm2.Fitted)*10;
nn = numel(pThr);
sens = zeros(1,nn);
spec = zeros(1,nn);

trueLab = DATA.ASPRESENCE_T72
for i=1:nn
    Ipred = lm2.Fitted>pThr(i);
    % P(ASPRESENCE|testPredictsPresence)
    sens(i) = sum(and(trueLab,Ipred))/sum(trueLab);
    spec(i) = sum(and(~trueLab,~Ipred))/sum(~trueLab);
    pause(.0)
end
clf
plot(spec,sens)





%% %% %% %% %% %%

[lm2,dev,stat] = mnrfit(X,Y);
xx = -[ones(height(DATA),1),X]*lm2;
p = 1./(1 + exp(-xx));
% compute probability of pathology given test is positive:


pThr = (0:.00001:0.01);
nn = numel(pThr);
sens = zeros(1,nn);
spec = zeros(1,nn);
for i=1:nn
    Ipred = p>pThr(i);
    % P(ASPRESENCE|testPredictsPresence)
    posPredVal =    sum(and(DATA.ARPRESENCE_T72,Ipred))/sum(Ipred);
    sens(i) = sum(and(DATA.ARPRESENCE_T72,Ipred))/sum(DATA.ARPRESENCE_T72);
    spec(i) = sum(and(~DATA.ARPRESENCE_T72,~Ipred))/sum(~DATA.ARPRESENCE_T72);
end


clf
plot(spec,sens)
stat.p