% In this script I wish to get a sense of what the noise looks like. Are
% they harder to segment?

% rename data frame for convenience:
data = HSdataTrain;
noiseData{1} = data.MURMUR_1NOISE_REF_T72;
noiseData{2} = data.MURMUR_2NOISE_REF_T72;
noiseData{3} = data.MURMUR_3NOISE_REF_T72;
noiseData{4} = data.MURMUR_4NOISE_REF_T72;



n = height(data);
clearvars N P I J noiseLvl

% in this script, I define a noisy recording as the observers agreeing on
% noise. 

% How many of the sound files are labelled as noisy?
I.noise1 = data.MURMUR_1NOISE_REF_T72; J.noise1 = find(I.noise1);
I.noise2 = data.MURMUR_2NOISE_REF_T72; J.noise2 = find(I.noise2);
I.noise3 = data.MURMUR_3NOISE_REF_T72; J.noise3 = find(I.noise3);
I.noise4 = data.MURMUR_4NOISE_REF_T72; J.noise4 = find(I.noise4);

noiseLvl{1} = data.MURMUR_1NOISE_REF_T72(J.noise1);
noiseLvl{2} = data.MURMUR_2NOISE_REF_T72(J.noise2);
noiseLvl{3} = data.MURMUR_3NOISE_REF_T72(J.noise3);
noiseLvl{4} = data.MURMUR_4NOISE_REF_T72(J.noise4);



% agreed noise in atleast one location
I.noiseAny = I.noise1+I.noise2+I.noise3+I.noise4>0;
J.noiseAny = find(I.noiseAny);

N.noiseAny = sum(I.noise1);
N.noisyTotal = sum(I.noise1)+sum(I.noise2)+sum(I.noise3)+sum(I.noise4);
P.noiseAny = N.noiseAny/n;
plotIt=false;
% we see that about 6%, roughly  1/17 of the recordings, are labelled as
% noisy. Not an insignificant amount. Lets plot the noisy recordings.
if plotIt
for i=1:N.noiseAny
    ii = J.noiseAny(i);
    id = data.UNIKT_LOPENR(ii);
    clf
    for aa=1:4
        x  = wav2TS(id,aa);
        x  = downsample(x,30);  
        subplot(2,2,aa)
            getScaleogram(x,1,true);
            noiseLvl = noiseData{aa}(ii);
            title(sprintf('noise=%.1g',noiseLvl))
    end
%     saveas(gcf,sprintf('noisyRec%g_id%.0f.png',i,id));
    pause(0)
end
end
% conclusion: The recordings that are labelled as noisy agree prety well
% with visual perception of noisiness judging from the scaleograms.
% However, there are some recordings labelled as noisy that do not look
% that bad. Undoubtedly, some of the recordings labelled noisy will be
% segmentable, and some of the labelled not noisy will be unsegmentable.

% the next question is; does excluding the noisy labels bias our results?
% In other words, is noisyness itself a feature that could be indicative of
% pathology?


I.AR = data.ARGRADE_T72>=2;
I.MR = data.MRGRADE_T72>=2;
I.AS = data.ASGRADE_T72>=1;
I.MS = data.MSGRADE_T72>=1;

P.AR = sum(I.AR)/n;
P.MR = sum(I.MR)/n;
P.AS = sum(I.AS)/n;
P.MS = sum(I.MS)/n;

P.noiseGivenAR     = condProb(I.noiseAny,I.AR);
P.noiseGivenMR     = condProb(I.noiseAny,I.MR);
P.noiseGivenAS     = condProb(I.noiseAny,I.AS);
P.noiseGivenMS     = condProb(I.noiseAny,I.MS);
P.noiseGivenRegurg = condProb(I.noiseAny,or(I.AR,I.MR));
P.ARgivenNoise     = condProb(I.AR, I.noiseAny);
P.MRgivenNoise     = condProb(I.MR, I.noiseAny);
% the prevalence of noise is higher for the AR population; it goes from 5.8%
% to 12%. Similar numbers hold for MR. Given AS, the prevalence of noise
% increases to 7.7%. In terms of changes to risk due to presence of noise,
% there does not seem to be any significant differences. So, the conclusion
% is that we can disregard noisy data without biasing our results; noise is
% not strongly correlated with any of the illnesses. This is the same
% conclusion that was found in the article we discussed together.


% Next quesiton: is the segmentation algorithm capable of handling noisy
% data?








