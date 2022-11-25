modelData = HSdata;

ASgrade = modelData.ASGRADE_T72;
MG      = [modelData.Murmur_1_grade_ref_ny_T72,modelData.Murmur_2_grade_ref_ny_T72,...
           modelData.Murmur_3_grade_ref_ny_T72,modelData.Murmur_4_grade_ref_ny_T72];
maxMG   = modelData.maxMeanMurGrade;
strokeVol = modelData.LVSTROKEVOL_T72;
ejectionFrac = modelData.LVEFBIPLANE_T72;
plot(ASgrade,maxMG + normrnd(0,1,[1,10]),'*')

%%
close all
figure
plot(maxMG./(100+strokeVol) + norm(1),ASgrade)
figure
boxplot(maxMG,ASgrade)
xlabel 'AS grade'
ylabel 'murmur grade'

%%
figure
plot(strokeVol)
figure
plot(ejectionFrac)
%%

clf
subplot(2,1,1)
boxplot(maxMG,ASgrade)
title(sprintf('pos. %g',i))
subplot(2,1,1)
boxplot(maxMG./strokeVol,ASgrade)
title(sprintf('pos. %g',i))




