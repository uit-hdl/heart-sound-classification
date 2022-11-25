atleastOneDiaMur = [HSdata.MURMUR_1DIA_REF_T72,HSdata.MURMUR_2DIA_REF_T72,...
                    HSdata.MURMUR_3DIA_REF_T72,HSdata.MURMUR_4DIA_REF_T72];
Idiastolic = sum(atleastOneDiaMur,2)>0;
Isystolic = and(HSdata.maxMeanMurGrade>0,~Idiastolic);
ARsig = HSdata.ARgrade>=3;
ARmur = and(ARsig,HSdata.maxMeanMurGrade>0);
IsystolicARmurmurs = and(ARmur,~Idiastolic);
ASorMRsig = or(HSdata.ASgrade>0,HSdata.MRgrade>=2);
systolicARmur_notMRorAS = and(IsystolicARmurmurs,~ASorMRsig);
sum(and(ARmur,~ASorMRsig))/sum(ARsig)

sum(and(~ASorMRsig,Isystolic))/sum(~ASorMRsig)
sum(and(~ASorMRsig,Isystolic))/height(HSdata)