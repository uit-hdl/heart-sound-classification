% some summary statistics
% mean = 27.1777
 %#ok<*NOPTS>
% it is annoying to write HSdataTrain all the time...
data = HSdataTrain;
data = data(~isnan(data.BMI_T7),:);
n = height(data);


% define some logical index vectors
clearvars I N P
thr.fat = 26;
I.fat  = data.BMI_T7>thr.fat;
I.slim = ~I.fat;
I.ARsick = data.ARGRADE_T72>=3;
I.MRsick = data.MRGRADE_T72>=3;
I.ASsick = data.ASGRADE_T72>=1;
I.MSsick = data.MSGRADE_T72>=1;
I.sick = (I.ARsick+I.MRsick+I.ASsick+I.MSsick)>0;

N.fat  = sum(I.fat);
N.slim = sum(I.slim);
N.sick = sum(I.sick); 
N.healthy = sum(~I.sick);
N.fatAndSick = sum(and(I.fat,I.sick));
meanBMI = mean(data.BMI_T7);

P.fat = Nfat/n;
P.slim = N.slim/n;
P.sick = Nsick/n;

% conditional probabilities
P.fatGivenSick = sum(and(I.fat,I.sick))/N.sick
P.slimGivenSick = sum(and(I.slim,I.sick))/N.sick
P.sickGivenFat = sum(and(I.sick,I.fat))/N.fat
P.sickGivenSlim = sum(and(I.sick,I.slim))/N.slim

% interestingly, the chance of being heart sick increases if you are slim.
% More accurately; the proportion of slim people who are sick is higher
% than the proportion of fat people who are sick. This might be because of
% the fact that old people tend to be underweight, and this acts as a
% confounding factor.

%%% Does BMI correlate negatively with age?

correl = corrcoef([data.BMI_T7,data.AGE_T7]);
% there is a very weak correlation of -0.0066. between age and BMI. At
% first glance, this does not seem to represent the difference very well.

correl = corrcoef([data.BMI_T7,data.MRGRADE_T72])
% interestingly there is a significant negative correlation of -0.1538
% between BMI and grade of MR. This agrees with the above observation that
% slim people have a higher proportion of heart sick. What is going on?

% One possibility is that the majority of the heartsick are very old, and
% that this population is quite small, so they are not well represented in
% the summary statistics. Lets have a look:

mean(data.AGE_T7(I.sick))
% the mean age of the heart sick is 70 years old!
thr.old = 81;
I.old = data.AGE_T7>thr.old;
I.young = ~I.old;
P.old = sum(I.old)/n;
P.slimGivenOld = sum(and(I.slim,I.old))/sum(I.old)
P.fatGivenOld = sum(and(I.fat,I.old))/sum(I.old)
P.oldGivenSlim = sum(and(I.old,I.slim))/sum(I.slim)
P.oldGivenFat = sum(and(I.old,I.fat))/sum(I.fat)
% A large difference. Of the people who are old, 75% of them are slim, and
% 25% are fat. But is this just reflecting the fact that there are more
% slim poeple by the way we have defined fat and slim?

mean(data.BMI_T7(I.old))
mean(data.BMI_T7(I.young))
% we see that the average BMI is lower for the old than for the you when
% the threshold is set to 75. The difference seems to skew more in the
% direction of the old being thinner as the threshold for "old" increases.

% Is the conclusion sensitive to the definition of old?
pslim = zeros(1,36);
pfat  = zeros(1,36);
psick = zeros(1,36);
for i=45:80
    I.old = data.AGE_T7>i;
    pslim(i-44) = sum(and(I.slim,I.old))/sum(I.old);
    pfat(i-44)  = sum(and(I.fat,I.old))/sum(I.old);
    psick(i-44) = sum(and(I.sick,I.old))/sum(I.old);
end
clf
plot(45:80,pslim,'b')
hold on
plot(45:80,pfat,'r')
plot(45:80,psick,'k')
title 'P(slim|age>x)'
ylim([0,1])
legend({'slim','fat'})

% We see that the chance of being sick rises quickly with age, and that as
% people pass the age of about 68, the chance of being slim starts rising.
% This may explain why the risk of being sick is slightly higher for the
% slim population. 

% What proportion of the sick are above the threshold for old age?
P.oldGivenSick = sum(and(I.sick,I.old))/sum(I.sick)
% Of the ones that are sick, 12.7% are old, yet the proportion of old
% people is only 3.6%. Clearly this age group is greatly overrepresented
% amongst the sick population, and since this group is also below average
% in weight, it makes sense that this would result in a low weight
% appearing to increase the risk of being sick.

% conclusions: using a threshold of 26 for being fat, and 65 for being old,
% we get that the chance of being sick given fat is 18.2%, the chance of
% being sick given slim is 20%, and the chance of being 

%%% ARE FAT PEOPLE HARDER TO AUSCULTATE?
% It remains to see if fat people are harder to auscultate due to increased
% levels of subcutaneous fat. The hypothesis is that
% P(murgrade>thr.mur|fat) should be lower than P(murgrade>thr.mur|slim),
% i.e. it should be harder to detect.

thr.mur = 2.0; % threshold for audible murmur.

I.mur = data.maxMeanMurGrade>thr.mur;
I.fatAndSick = and(I.fat,I.sick);
I.slimAndSick = and(I.slim,I.sick);
P.murGivenFatAndSick = sum(and(I.mur,I.fatAndSick))/sum(I.fatAndSick)
P.murGivenSlimAndSick = sum(and(I.mur,I.slimAndSick))/sum(I.slimAndSick)

ci.murGivenFatAndSick  = ciPest(I.mur,I.slimAndSick)
ci.murGivenSlimAndSick = ciPest(I.mur,I.fatAndSick)
% we get that the chance to detect sickness for fat people is 16.1% and the
% chance to detect sickness for slim people is 19.6%. It is not entirely
% clear if this is significant. the confidence intervals for these
% estimates are [0.1320 0.2599] and [0.1063 0.2155] respectively, not a
% significant difference.

% Conclusion: There is not evidence that BMI significantly alters the
% probability of detecting valvular heart disease through auscultation.


% Finally, lets see if there is a significant difference between the
% conditional risks of being sick considering BMI:

ci.sickGivenSlim = ciPest(I.sick,I.slim)
ci.sickGivenFat = ciPest(I.sick,I.fat)

% sickGivenSlim: [0.1655 0.2310], sickGivenFat: [0.1629 0.2083]
% Not very significant

% How sensitive is this conclusion with regards to the chosen threshold for
% fat?
psick = zeros(20,32);
for i=20:32
    Ifat = data.BMI_T7>i;
    Islim = ~Ifat;
    pRiskFat(i-19) = sum(and(I.sick,Ifat))/sum(Ifat);
    pRiskSlim(i-19) = sum(and(I.sick,Islim))/sum(Islim);
end
clf
plot(20:32,pRiskSlim,'b')
hold on
plot(20:32,pRiskFat,'r')
title 'conditional risks for slim and fat population respectively as function of threshold'
ylim([0,.5])
legend({'slim','fat'})

% the conclusions do not seem to be very sensitive to choice of obesity
% threshold.




