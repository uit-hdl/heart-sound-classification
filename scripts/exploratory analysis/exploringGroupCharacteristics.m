% First I want to explore the subgroup consisting of those with AR>=3. What
% variables characterize them?

Iall = 1:height(HSdata0);
Ind  = stats.I.AR.GEQ.grade3.any;
data = HSdata0(Ind,:);

%% what is the gender distribution for sick cohort?
mean(data.SEX_T7==1) % 55.63% male
mean(data.SEX_T7==0) % 44.37% female
%% what is the gender distribution for sick cohort?
mean(data.AGE_T7)    % 71.3113 
mean(HSdata0.AGE_T7) % 64.0005
% 71.3113 vs 64.0005
%% what is the 

%%


histogram(HSdata0(Ind,:).SEX_T7,'normalization','probability')
subplot(211)
histogram(HSdata0(Iall,:).DYSPNEA_FAST_UPHILL_T7,'normalization','probability')
subplot(212)
histogram(HSdata0(Ind,:).DYSPNEA_FAST_UPHILL_T7,'normalization','probability')

clf
subplot(211)
    histogram(HSdata0(Iall,:).PO2_T72,'normalization','probability')
    xlim([90,102])
subplot(212)
    histogram(HSdata0(Ind,:).PO2_T72,'normalization','probability')
    xlim([90,102])
    
    subplot(211)
        histogram(HSdata0(Iall,:).CHEST_PAIN_TIME_T7,'normalization','probability')
    subplot(212)
        histogram(HSdata0(Ind,:).CHEST_PAIN_TIME_T7,'normalization','probability')

subplot(211)
    histogram(HSdata0(Iall,:).PULSESPIRO_T72	,'normalization','probability')
    xlim([30,110])
subplot(212)
    histogram(HSdata0(Ind,:).PULSESPIRO_T72,'normalization','probability')
    xlim([30,110])



% DYSPNEA_CALMLY_FLAT not very informative
% DYSPNEA_FAST_UPHILL very informative, about half have it
% CHEST_PAIN does not seem very significant
