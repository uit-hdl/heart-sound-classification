murDataInfo.allMurData.col = (41:146);

% observer 1
murDataInfo.classifObs1.col = (41:56);
murDataInfo.classifObs1.description = sprintf('variables describe obs.1 classification of the heart sounds for each of the 4 different auscultation locations. \n Gives information about: \nAuscultation Location (1,2,3,4) \n Heart Cycle Location (systole, diastole) \n Sound type (normal, murmur, noise)');
murDataInfo.classifObs1.nameFormat = 'Murmur_{1/2/3/4}{sys/dia}AD_T72';


murDataInfo.gradeObs1.col = (57:60);
murDataInfo.gradeObs1.descripton = sprintf('grading of the murmurs by obs. 1 from 1-6 for each Ausultation Location.');
murDataInfo.gradeObs1.nameFormat = 'Murmur_{1/2/3/4}gradeNRAD_T72';

murDataInfo.systMurAndS2Obs1.col = (61:72);
murDataInfo.systMurAndS2Obs1.descripton = sprintf('Indicator variable that assumes presence of murmur. Desribes whether or \nnot, according to obs. 1, S2 and systolic murmur is audible, and whether or not \nS2 is seperable from systolic murmur. Information is given for each \nof the 4 Auscultation Locations.');
murDataInfo.systMurAndS2Obs1.nameFormat = 'Murmur_2tone_{1/2/3/4}{na/usep/sep}AD_T72';

% observer 2
murDataInfo.classifObs2.col = (73:88);
murDataInfo.classifObs2.description = sprintf('variables describe obs. 2 classification of the heart sounds for each of the 4 different auscultation locations. \n Gives information about: \nAuscultation Location (1,2,3,4) \n Heart Cycle Location (systole, diastole) \n Sound type (normal, murmur, noise)');
murDataInfo.classifObs2.nameFormat = 'Murmur_{1/2/3/4}{sys/dia}SA_T72';

murDataInfo.gradeObs2.col = (89:92);
murDataInfo.gradeObs2.descripton = sprintf('grading of the murmurs by obs. 2 from 1-6 for each Ausultation Location.');
murDataInfo.gradeObs2.nameFormat = 'Murmur_{1/2/3/4}gradeNRSA_T72';

murDataInfo.systMurAndS2Obs2.col = (93:103);
murDataInfo.systMurAndS2Obs2.descripton = sprintf('Indicator variable that assumes presence of murmur. Desribes whether or \nnot, according to obs. 2, S2 and systolic murmur is audible, and whether or not \nS2 is seperable from systolic murmur. Information is given for each \nof the 4 Auscultation Locations.');
murDataInfo.systMurAndS2Obs2.nameFormat = 'Murmur_2tone_{1/2/3/4}{na/usep/sep}SA_T72';


% Murmur Agreement Variables
murDataInfo.soundAgree.col = (104:119);
murDataInfo.soundAgree.description = sprintf('Desribes the agreement between obs. 1 and 2 on the type of sound heard \n(normal, murmur, noise) for each of the 4 Auscultation Locations. 1 indicates \nagreement and 0 disagreement.');
murDataInfo.soundAgree.nameFormat = 'Murmur_{1/2/3/4}{norm/sys/dia/noise}_ref_T72}';

% Murmur carotis
murDataInfo.carotidMurmur.col = (120:121);
murDataInfo.carotidMurmur.description = sprintf('escribes whether or not a continuous murmur is heard, which indicates \ncarotid artery insufficiency.');
murDataInfo.carotidMurmur.nameFormat = 'Murmur_{carl/carr}_ref_T72';

% Auidibility of S1 acc. to obs. 1
murDataInfo.auidibilityS1Obs1.col = (122:129);
murDataInfo.auidibilityS1Obs1.description = sprintf('Describes how audible S1 is according to obs. 1, rated either faint or \ninaudible.');
murDataInfo.auidibilityS1Obs1.nameFormat = 'First_sound_{1/2/3/4}{faint/inaudible}AD_T72';

% Auidibility of S1 acc. to obs. 2
murDataInfo.auidibilityS1Obs2.col = (130:137);
murDataInfo.auidibilityS1Obs2.description = sprintf('Describes how audible S1 is according to obs. 2, rated either faint or \ninaudible for each auscultation location.');
murDataInfo.auidibilityS1Obs2.nameFormat = 'First_sound_{1/2/3/4}{faint/inaudible}SA_T72';


% Presence of murmurs in atleast one auscultation location
murDataInfo.atLeastOneSystMurm.col = (138:140);
murDataInfo.atLeastOneSystMurm.description = sprintf('Could a systolic murmur be heard in atleast one auscultation location? \n Information from each observer, as well as agreement score.');
murDataInfo.atLeastOneSystMurm.nameFormat  = 'Murmur_sys_1234_{AD/SA/ref}_T72';

murDataInfo.NlocSoundAgree.col = (141:144);
murDataInfo.NlocSoundAgree.description = sprintf('In how many locations was there agreement on type of sound? Takes \nvalues from 0 to 4. ');
murDataInfo.NlocSoundAgree.nameFormat  = '{Normal/Syst/Dia/Noise}_number_1234_ref_T72';

murDataInfo.allNormalOrNoise.col = (145:145);
murDataInfo.allNormalOrNoise.description = sprintf('1 if only normal sounds or noise could be heard.');
murDataInfo.allNormalOrNoise.nameFormat = 'Normal_total_T72';

murDataInfo.isClassifiable.col = (146:146);
murDataInfo.isClassifiable.description = sprintf('Is the sound classifiable? Sound is defined as classifiable if either holds: \n* sound is classified as normal/noise in <=2 auscultation locations\n* systolic or diastolic mumrmur has been agreed upon for atleast one ausultation location.');
murDataInfo.isClassifiable.nameFormat = 'Classifyable_set_T72';

%% Generate missing data; mean murmur grade for each ausc. location
T = 1/2*(murData{:,murDataInfo.gradeObs1.col} + murData{:,murDataInfo.gradeObs2.col});

murData.Murmur_1_grade_ref_ny_T72 = T(:,1);
murData.Murmur_2_grade_ref_ny_T72 = T(:,2);
murData.Murmur_3_grade_ref_ny_T72 = T(:,3);
murData.Murmur_4_grade_ref_ny_T72 = T(:,4);
clear T

% Add descriptions of mean murmur grade data column
murDataInfo.meanMurGrade.col = (258:261);
murDataInfo.allMurData.col = [murDataInfo.allMurData.col,(258:261)];
murDataInfo.meanMurGrade.description = sprintf('Mean murmur grade (1-6) for each of the 4 auscultation locations.');
murDataInfo.meanMurGrade.nameFormat = 'Murmur_{1/2/3/4}_grade_ref_ny_T72';
murDataInfo.meanMurGrade.newData = true;

%% Generate missing data; Distinct murmur any location
murDa = murDataInfo.meanMurGrade.col;
murData.Murmur_1234_distinct_T72 = sum(murData{:,murDa}>=2.5,2)>0;

murDataInfo.murDistinct.col = (262:262);
murDataInfo.allMurData.col(end+1) = 262;
murDataInfo.murDistinct.description = sprintf('indicates whether or not a distinct murmur could be heard in atleast one \nauscultation location. A murmur is defined as distinct when its mean grade is \natleast 2.5.');
murDataInfo.murDistinct.nameFormat = 'Murmur_1234_distinct_T72';
murDataInfo.murDistinct.newData = true;
%% Generate missing data; Distinct or faint murmur any location
murDa = murDataInfo.meanMurGrade.col;
murData.Murmur_1234_faint_T72 = sum(...
    and( murData{:,murDa}>=0.5, murData{:,murDa}<2.5) ,2)...
                                        >0;

murDataInfo.murFaint.col = (263:263);
murDataInfo.allMurData.col(end+1) = 263;
murDataInfo.murFaint.description = sprintf('indicates whether or not a faint murmur could be heard in atleast one \nauscultation location. A murmur is defined as faint when its mean grade is \nbetween 0.5 and (less than) 2.5.');
murDataInfo.murFaint.nameFormat = 'Murmur_1234_faint_T72';
murDataInfo.murFaint.newData = true;

%% Generate missing data; Distinct OR faint murmur any location
murData.Murmur_1234_faint_distinct_T72 = ...
    or(murData.Murmur_1234_faint_T72,murData.Murmur_1234_distinct_T72);

murDataInfo.murFaintDistinct.col = (264:264);
murDataInfo.allMurData.col(end+1) = 264;
murDataInfo.murFaintDistinct.description = sprintf('indicates whether or not a faint or distinct murmur could be heard in atleast one \nauscultation location. A murmur is defined as distinct or faint when its mean grade is \natleast 0.5.');
murDataInfo.murFaintDistinct.nameFormat = 'Murmur_1234_faint_distinct_T72';
murDataInfo.murFaintDistinct.newData = true;
%% Generate missing data; Distinct OR faint murmur any location
murDa = murDataInfo.meanMurGrade.col;
m = max(murData{:,murDa},[],2);
murData.maxMurAuscLoc = cell(size(murData,1),1);
for i=1:size(murData,1)
    murData.maxMurAuscLoc{i} = find(murData{i,murDa}==m(i));
end

murDataInfo.maxMurAuscLoc.col = (265:265);
murDataInfo.allMurData.col(end+1) = 265;
murDataInfo.maxMurAuscLoc.description = sprintf('Gives location of the auscultation location where mean murmur grade is highest. \nif there is a tie the locations are given in vector form.');
murDataInfo.maxMurAuscLoc.nameFormat = 'maxMurAuscLoc';
murDataInfo.maxMurAuscLoc.newData = true;
%% Generate Data: maximum mean murmur grade
murDa = murDataInfo.meanMurGrade.col;
murData.maxMeanMurGrade = max(murData{:,murDa},[],2);

murDataInfo.maxMeanMurGrade.col = (266:266);
murDataInfo.allMurData.col(end+1) = 266;
murDataInfo.maxMeanMurGrade.description = sprintf('Gives the maximum mean murmur grade.');
murDataInfo.maxMeanMurGrade.nameFormat = 'maxMeanMurGrade';
murDataInfo.maxMeanMurGrade.newData = true;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Create tables that are handy for getting visual overwiev of murmur data
murDataInfo.SummaryTableObs1 = cell(size(murData,1),1);
murDataInfo.SummaryTableObs2 = cell(size(murData,1),1);
murDataInfo.SummaryTableObs12 = cell(size(murData,1),1);

for k=1:size(murData,1)
    temp = murData{k,murDataInfo.classifObs1.col}';

    normal    = temp(1:4:end);
    systolic  = temp(2:4:end);
    diastolic = temp(3:4:end);
    noise     = temp(4:4:end);
    murmurGrade = murData{k,murDataInfo.gradeObs1.col}';

    summaryTable = table(["AusLoc 1";"AusLoc 2";"AusLoc 3";"AusLoc 4"],...
                    normal, systolic, diastolic, noise, murmurGrade);
    murDataInfo.SummaryTableObs1{k} = summaryTable;
end

for k=1:size(murData,1)
    temp = murData{k,murDataInfo.classifObs2.col}';

    normal    = temp(1:4:end);
    systolic  = temp(2:4:end);
    diastolic = temp(3:4:end);
    noise     = temp(4:4:end);
    murmurGrade = murData{k,murDataInfo.gradeObs2.col}';

    summaryTable = table(["AusLoc 1";"AusLoc 2";"AusLoc 3";"AusLoc 4"],...
                    normal, systolic, diastolic, noise, murmurGrade);
    murDataInfo.SummaryTableObs2{k} = summaryTable;
end

for k=1:size(murData,1)
    temp1 = murData{k,murDataInfo.classifObs1.col}';
    temp2 = murData{k,murDataInfo.classifObs2.col}';

    normal    = [temp1(1:4:end),temp2(1:4:end)];
    systolic  = [temp1(2:4:end),temp2(2:4:end)];
    diastolic = [temp1(3:4:end),temp2(3:4:end)];
    noise     = [temp1(4:4:end),temp2(4:4:end)];
    murmurGrade = [murData{k,murDataInfo.gradeObs1.col}',...
                    murData{k,murDataInfo.gradeObs2.col}'];
    murmurGradeMean = murData{k,murDataInfo.meanMurGrade.col}';
    
    summaryTable = table(["AusLoc 1";"AusLoc 2";"AusLoc 3";"AusLoc 4"],...
                    normal, systolic, diastolic, noise, murmurGrade, murmurGradeMean);
    murDataInfo.SummaryTableObs12{k} = summaryTable;
end


