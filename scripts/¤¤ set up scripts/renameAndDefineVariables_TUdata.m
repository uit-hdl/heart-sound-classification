%% simplify some notation and add variables:
HSdata.id = HSdata.UNIKT_LOPENR;
HSdata.murGrade1 = HSdata.Murmur_1_grade_ref_ny_T72;
HSdata.murGrade2 = HSdata.Murmur_2_grade_ref_ny_T72;
HSdata.murGrade3 = HSdata.Murmur_3_grade_ref_ny_T72;
HSdata.murGrade4 = HSdata.Murmur_4_grade_ref_ny_T72;

HSdata.noise1 = HSdata.MURMUR_1NOISE_REF_T72==1;
HSdata.noise2 = HSdata.MURMUR_2NOISE_REF_T72==1;
HSdata.noise3 = HSdata.MURMUR_3NOISE_REF_T72==1;
HSdata.noise4 = HSdata.MURMUR_4NOISE_REF_T72==1;

HSdata = renamevars(HSdata,'ARGRADE_T72','ARgrade');
HSdata = renamevars(HSdata,'MRGRADE_T72','MRgrade');
HSdata = renamevars(HSdata,'ASGRADE_T72','ASgrade');
HSdata = renamevars(HSdata,'MSGRADE_T72','MSgrade');

HSdata = renamevars(HSdata,'DIABETES_T7','diabetes');
HSdata = renamevars(HSdata,'DYSPNEA_CALMLY_FLAT_T7','dyspneaCalmlyFlat');
HSdata = renamevars(HSdata,'DYSPNEA_FAST_UPHILL_T7','dyspneaFastUpphill');
HSdata = renamevars(HSdata,'DYSPNOE_REST_T7','dyspneaRest');
HSdata = renamevars(HSdata,'ANGINA_T7','angina');
HSdata = renamevars(HSdata,'PO2_T72','po2');
HSdata = renamevars(HSdata,'PULSESPIRO_T72','pulseSpiro');
HSdata = renamevars(HSdata,'CHEST_PAIN_NORMAL_T7','chestPainNormal');
HSdata = renamevars(HSdata,'CHEST_PAIN_FAST_T7','chestPainFast');
HSdata = renamevars(HSdata,'CHEST_PAIN_ACTION_T7','chestPainAction');
HSdata = renamevars(HSdata,'AGE_T7','age');
HSdata = renamevars(HSdata,'BMI_T7','bmi');
HSdata = renamevars(HSdata,'SEX_T7','sex');
HSdata = renamevars(HSdata,'AVMEANPG_T72','avmeanpg');
HSdata = renamevars(HSdata,'AVAVMAX_T72','avarea');
HSdata = renamevars(HSdata,'SMOKE_DAILY_Q2_T7','smoke');

HSdata.dyspneaRestOrFlat = myor(HSdata.dyspneaRest, HSdata.dyspneaCalmlyFlat,"and");
HSdata.anginaOrDyspnea   = myor(HSdata.angina, HSdata.dyspneaRestOrFlat,"and");
HSdata.chestPain = myor(HSdata.chestPainNormal, HSdata.chestPainFast,"and");
HSdata.highBP    = (HSdata.HIGH_BLOOD_PRESSURE_T7>0) + nanVec(HSdata.HIGH_BLOOD_PRESSURE_T7);

HSdata.ARsigSympt = and(HSdata.ARgrade>=3,HSdata.anginaOrDyspnea>0);
HSdata.MRsigSympt = and(HSdata.MRgrade>=3,HSdata.anginaOrDyspnea>0);

HSdata.smokeCurrent = (HSdata.smoke==1) + nanVec(HSdata.smoke); % combine previous and current smokers
HSdata.smoke    = or(HSdata.smoke==1,HSdata.smoke==2) + nanVec(HSdata.smoke); % combine previous and current smokers
HSdata.angina   = (HSdata.angina>0) + nanVec(HSdata.angina); % combine those that have or have had angina petctoris
HSdata.diabetes = (HSdata.diabetes>0) + nanVec(HSdata.diabetes); % combine those who have or have had diabetes

% cases where all positions had noisy audio:
HSdata.noiseOnly = sum([HSdata.noise1,HSdata.noise2,HSdata.noise3,HSdata.noise4],2)==4;

% diastolic murmur in atleast 1 position:
HSdata.diastMurPresence = myor({HSdata.MURMUR_1DIA_REF_T72,...
                                HSdata.MURMUR_2DIA_REF_T72,...
                                HSdata.MURMUR_3DIA_REF_T72,...
                                HSdata.MURMUR_4DIA_REF_T72});
                           

HSdata.murGradeSum = sum([HSdata.murGrade1,HSdata.murGrade2,HSdata.murGrade3,HSdata.murGrade4],2);
HSdata.murGradeMax = max([HSdata.murGrade1,HSdata.murGrade2,HSdata.murGrade3,HSdata.murGrade4],[],2);

%% Redefining AS from mean-pressure-gradient:
% note: if avmeanpg was not available, then the original AS grade was used
% as as the grade for the participant. There were 23 such cases.
HSdata.ASgrade0 = HSdata.ASgrade;

Iavmpg_none = HSdata.avmeanpg<15;
Iavmpg_mild = and(15<=HSdata.avmeanpg,HSdata.avmeanpg<20);
Iavmpg_moderate = and(20<=HSdata.avmeanpg,HSdata.avmeanpg<40);
Iavmpg_severe = and(40<=HSdata.avmeanpg,HSdata.avmeanpg<inf);

HSdata.ASgrade(Iavmpg_none) = 0;
HSdata.ASgrade(Iavmpg_mild) = 1;
HSdata.ASgrade(Iavmpg_moderate) = 2;
HSdata.ASgrade(Iavmpg_severe)   = 3;

%% add column with indicator variable for significant VHD (def. by grade)
T = [HSdata.ARgrade>=3,HSdata.MRgrade>=3,HSdata.ASgrade>=1,HSdata.MSgrade>=1];
T = sum(T,2)>0;
HSdata.sigVHD31 = T;

T = [HSdata.ARgrade>=4,HSdata.MRgrade>=4,HSdata.ASgrade>=1,HSdata.MSgrade>=1];
T = sum(T,2)>0;
HSdata.sigVHD41 = T;
%% add column with indicator variable for significant VHD (def. by grade and symptoms)
T = [HSdata.ARsigSympt>0,HSdata.MRsigSympt>0,HSdata.ASgrade>=1,HSdata.MSgrade>=1];
T = sum(T,2)>0;
HSdata.sigSymptVHD = T;
%% add column with indicator variable for significant assymptomatic VHD:
HSdata.sigAsymptVHD = myor({and(HSdata.ARgrade>2, HSdata.anginaOrDyspnea==0),...
                         and(HSdata.MRgrade>2, HSdata.anginaOrDyspnea==0),...
                         and(HSdata.ASgrade>0, HSdata.anginaOrDyspnea==0),...
                         and(HSdata.MSgrade>0, HSdata.anginaOrDyspnea==0)});

%% get rows with atleast one usable (non-noisy) recording:
I = unionIterated({HSdata.noise1==0,HSdata.noise2==0,...
                    HSdata.noise3==0,HSdata.noise4==0},"logical");
HSdata.atleastOneClean = I;