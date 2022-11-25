
echoDataInfo.allEchoData.col = (147:257);

echoDataInfo.insuffAndSclerosis.col = (153:156);
echoDataInfo.insuffAndSclerosis.description = sprintf('indicators variables for AS, AI, MS, MI.');
echoDataInfo.insuffAndSclerosis.nameFormat = '{MI/AS/AI/MS}_T72';

echoDataInfo.stenosisAndRegurgPresence.col = [177,179,182,183,185];
echoDataInfo.stenosisAndRegurgPresence.description = sprintf('Indicates presence or absence of \n*Aortic Regurgitation\n*Mitral Regurgitation\n*Tricuspid Regurgitation\n*Mitral Stenosis\n*Atrial Stenosis');
echoDataInfo.stenosisAndRegurgPresence.nameFormat = '{AR/MR/TR/MS/AS}PRESENCE_T72';

echoDataInfo.stenosisAndRegurgGrade.col = [178,180,184,186];
echoDataInfo.stenosisAndRegurgGrade.description = sprintf('grades for \n*Aortic Regurgitation (trace, mild, moderate, severe)\n*Mitral Regurgitation (trace, mild, moderate, severe)\n*Mitral Stenosis (mild, average, severe)\n*Aortic Stenosis (mild, average, severe).');
echoDataInfo.stenosisAndRegurgGrade.nameFormat = '{AR/MR/MS/AS}GRADE_T72';

echoDataInfo.stenAndRegDecisionVars.col = [212,218];
echoDataInfo.stenAndRegDecisionVars.description = sprintf('Continuous variables used to decide grade of the 4 main types of VHD.\n AS -- aortic regurgitation flow pressure maximum gradient (*>30 --> AS)\n');

echoDataInfo.aorticValveData.col = 233:241;
echoDataInfo.aorticValveData.description = sprintf('Variables that describes the state of the Aortic Valve. Relevant for diagnosis of Aortic Stenosis.\nVariables included are\n*aortic Valve flow max pressure gradient\n*aortic valve flow mean pressure gradient\n*aortic valve flow maximum velocity\n*aortic valve flow mean velocity\n*aortic valve time-velocity\n*aortic valve area');




