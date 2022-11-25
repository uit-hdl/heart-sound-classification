% In this script I study the associations between murmur and VHD. 
clearvars I N

data = HSdataTrainAndVal;

% how many diastolic murmurs can be heard?
I.diast = data.MURMUR_4DIA_REF_T72; 
I.mur = data.Murmur_4_grade_ref_ny_T72; I.mur = I.mur>=1;
I.AI = data.AI_T72;
N.diast = sum(I.diast);
N.mur = sum(I.mur);
% only 7!
I.AR = data.ARGRADE_T72>=3; 
I.MR = data.MRGRADE_T72>=3; 
N.AR = sum(I.AR);
% 150 of degree moderate to severe. CleARly, the vast majority ARe not
% diast murs. So what is the mur that we ARe heARing I wonder?
% One possibility is that what we ARe heARing is AS, and that AR is just
% piggibacking on AS. How many of those with AR also have AS?

% How many of the diast sounds heARd at the mitral valve belong to AR?
N.ARandDiast = sum(and(I.diast,I.AR));
% of the 7 diast murs, 4 belonged to AR. Oddly, we heAR 5 diast
% murs in the aortic position, despite the fact that this is not where
% AR associated mur is supposed to be heARd well.
I.ARandMur = and(I.AR,I.mur);
I.ARandNoMur = and(I.AR,~I.mur);
I.ARandDiast = and(I.AR,I.diast);
N.ARandMur = sum(I.ARandMur);
% of the 150 cases of AR, 29 had murs, and 25 of these were not
% diast. How many ARe accounted for by the association with AS?
I.AS = data.ASGRADE_T72>0;
N.AS = sum(I.AS);
I.ASandAR = and(I.AS,I.AR);
N.ASandAR = sum(I.ASandAR);
% of 142 cASes of AR, 15 have AS. How many of the mur cASes of AR can be
% attributed to AS?
N.ASandARandMur = sum(and(I.AS,I.ARandMur));
% of the 29 murs, 12 were cASes of AS, meaning that is most likely what
% we heARd. That leaves us with 17 cASes of AR with mur that can not be
% attributed to AS. How many cASes ARe there of AR systolic murs that
% can not be attributed to AS?
I.ARandSysMurAndNotAS = and(not(I.AS),I.ARandMur).*not(I.diast);
N.ARandSysMurAndNotAS = sum(and(not(I.AS),I.ARandMur).*not(I.diast));
% Of these 17 murs that can not be attributed to AS, 12 were systolic.

% How many of the MS cases have diastolic murmurs?
I.MS = data.MSGRADE_T72>0; 
N.MS = sum(I.MS,'omitnan');
N.MSandDiast = sum(I.MS.*I.diast,'omitnan');
% surprisingly, none of the diastolic murmurss belonged to MS, and it does
% not matter which location we look at.
N.MSandMur = sum(I.mur.*I.MS,'omitnan');
% 11 of the 13 cases of MS had murs, and none of these were diastolic.

%%
labels.mur = true;
labels.vhd = true;
aa = 4;
k = k+1;
M = getScaleogramFromIndex(data,[],k,aa,true,nodes,20,labels);
% states = getExpandedReprOfStates(nodes.loc{1,aa},nodes.state{k,aa},n)
% there are indeed cases of AR where a clear systolic murmur can be
% detected which can not be attributed to AS. Hasse seems to be right

% what do the diastolic murmurs look like?
I.MR = data.MRGRADE_T72>0;



