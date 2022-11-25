% In this script, I explore to what degree we can expect AR to be predicted
% by murmurs. This is especially relevant now since I have just developed
% a murmur detection algorithm, and want to use it to see how good it is at
% detecting AR of different degrees. In order to get a sense of this, I
% need to have something to compare against as a base line. An interesting
% question is whether or not the algorithm is able to predict AR as good or
% better than the murmurs.

%#ok<*NOPTS>
clearvars I data

% decision variables:
aa = 1;
data = HSdataVal;
murthr = 2;
sickthr = 3;
VHD =  'ARGRADE_T72';

% define some useful variables:
murStr = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
I.mur  = data.(murStr)>=murthr;
I.sick = data.(VHD)>=sickthr;
N.mur = sum(I.mur);

% absolute probabilities
P.sick.null = prob(I.sick);
P.mur.null  = condProb(I.mur,[]); 
P.nomur.null = condProb(not(I.mur),[]);

% conditional probabilities
P.mur.sick = condProb(I.mur,I.sick);
P.sick.mur = condProb(I.mur,I.sick);

P.nomur.nosick = condProb(not(I.mur),not(I.sick));
P.nosick.nomur = condProb(not(I.sick),not(I.mur));

% characteristics of test that predicts AR if murmur grade is atleast 2,
% and no AR if not.
testchar.sens = P.mur.sick;
testchar.spec = P.nomur.nosick;
testchar.accu = P.mur.sick*prob(I.mur) + P.nomur.nosick*prob(~I.mur)

testSens =@(psi) psi + P.mur.sick*(1-psi);
testSpec =@(psi) (1-psi)*P.nomur.nosick;

%% plot the test:
psivec = -1:.01:1;

clf
plot(testSens(psivec),testSpec(psivec))
hold on
yline(1)
xline(0)
xline(testchar.sens)
yline(testchar.spec)

