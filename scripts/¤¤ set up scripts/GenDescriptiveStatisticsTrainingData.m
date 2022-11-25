data = HSdata0;
% how many have detactable murmur?
stats.I.audMur.any = data.maxMeanMurGrade>=2;
stats.N.any = height(data);
stats.N.audMur.any = sum(stats.I.audMur.any);
definitions.audMur = sprintf('The maximum of the mean murmur grades is greater than or equal to 0.5.\nNote that this is a very weak definition, and requires only one observer to hear murmur in location.');

% ROWS WITH DIASTOLIC MURMURS
IdiaMurObs1 = 0;
IdiaMurObs2 = 0;
for i=1:4
    str1 = sprintf('MURMUR_%gDIAAD_T72',i)
    str2 = sprintf('MURMUR_%gDIASA_T72',i)
    IdiaMurObs1 = IdiaMurObs1 + data.(str1);
    IdiaMurObs2 = IdiaMurObs2 + data.(str2);
end

% Diastolic murmurs
stats.genData.maxDiaMurGrade = max(IdiaMurObs1,IdiaMurObs2);
stats.I.diaMur.Obs1.any = IdiaMurObs1>0;
stats.I.diaMur.Obs2.any = IdiaMurObs2>0;
stats.I.diaMur.agree.any  = min(stats.I.diaMur.Obs1.any,stats.I.diaMur.Obs2.any);
stats.I.diaMur.obs1or2.any = max(stats.I.diaMur.Obs1.any,stats.I.diaMur.Obs2.any);
stats.N.diamur.agree.any = sum(stats.I.diaMur.agree.any);
stats.N.diamur.obs1or2.any = sum(stats.I.diaMur.obs1or2.any);


% AR,MR,AS,MS indicator vectors and numbers for different grades
VHDnames = {'AR','MR','AS','MS'};

for i=1:4
    VHDStr  = VHDnames{i};
    stats.I.(VHDnames{i}).presence.any = data.(strcat(VHDnames{i},'PRESENCE_T72'));
    stats.N.(VHDnames{i}).presence.any = sum(stats.I.(VHDnames{i}).presence.any,'omitnan');
    % cycle through grades
    for j=0:4
        if i>2 && j==4
            continue
        end

        gradeStr= sprintf('grade%g',j);
        stats.I.(VHDStr).GEQ.(gradeStr).any = (data.(strcat(VHDStr,'GRADE_T72'))>=j);
        stats.I.(VHDStr).EQ.(gradeStr).any  = (data.(strcat(VHDStr,'GRADE_T72'))==j);
        % how many with pathology have grade of so-and-so?
        stats.N.(VHDStr).GEQ.(gradeStr).any = sum(stats.I.(VHDStr).GEQ.(gradeStr).any);
        stats.N.(VHDStr).EQ.(gradeStr).any  = sum(stats.I.(VHDStr).EQ.(gradeStr).any);
        
        % how many with PATH of GRADE (relation to) x had AUDIBLE MURMUR?
        stats.N.(VHDStr).GEQ.(gradeStr).audMur = sum(and(stats.I.(VHDStr).GEQ.(gradeStr).any,...
                                                             stats.I.audMur.any));
        stats.N.(VHDStr).GEQ.(gradeStr).noAudMur = sum(and(stats.I.(VHDStr).GEQ.(gradeStr).any,...
                                                              ~stats.I.audMur.any));
                                                          
        
        stats.N.(VHDStr).EQ.(gradeStr).audMur = sum(and(stats.I.(VHDStr).EQ.(gradeStr).any,...
                                                             stats.I.audMur.any));
        stats.N.(VHDStr).EQ.(gradeStr).noAudMur = sum(and(stats.I.(VHDStr).EQ.(gradeStr).any,...
                                                           ~stats.I.audMur.any));
        
        % P(PATH of GRADE >= x | AUDIBLE MURMUR)?                                       
        stats.P.(VHDStr).GEQ.(gradeStr).audMur = stats.N.(VHDStr).GEQ.(gradeStr).audMur/...
                                                   stats.N.audMur.any ;
        % P(AUDIBLE MURMUR | PATH of GRADE >= x)?                               
        stats.P.audMur.(VHDStr).GEQ.(gradeStr) = stats.N.(VHDStr).GEQ.(gradeStr).audMur/...
                                                   stats.N.(VHDStr).GEQ.(gradeStr).any ;
        
        % P(PATH of GRADE == x | AUDIBLE MURMUR)?                                       
        stats.P.(VHDStr).EQ.(gradeStr).audMur = stats.N.(VHDStr).EQ.(gradeStr).audMur/...
                                                   stats.N.audMur.any ;
        % P(AUDIBLE MURMUR | PATH of GRADE == x)?                               
        stats.P.audMur.(VHDStr).EQ.(gradeStr) = stats.N.(VHDStr).EQ.(gradeStr).audMur/...
                                                   stats.N.(VHDStr).EQ.(gradeStr).any ; 
                                                         
    end
end

definitions.GEQ = sprintf('Greater than or Equal to.');

% clinically significant heart disease
stats.I.clinSignVHD.any = (stats.I.AR.GEQ.grade3.any + stats.I.MR.GEQ.grade3.any +...
            stats.I.AS.GEQ.grade1.any+stats.I.MS.GEQ.grade1.any)>=1;
stats.N.clinSignVHD.any = sum(stats.I.clinSignVHD.any);
stats.P.clinSignVHD.any = stats.N.clinSignVHD.any/stats.N.any;
definitions.clinSignVHD = sprintf('refers to subjects with \nAR>=3 \nMR>=3 \nAS>=1 \nMS>=1');

% innocent murmurs
stats.I.innocentMurm.any = and(~stats.I.clinSignVHD.any,stats.I.audMur.any);
stats.N.innocentMurm.any = sum(stats.I.innocentMurm.any);
%% CREATE DESCRIPTIVE TABLES
clear T
NAMES = {'AR','MR','AS','MS'};
% cells to be converted into tables:
T{1} = nan(6,4);
T{2} = nan(6,4);
T{3} = nan(6,4);
T{4} = nan(6,4);
T{5} = nan(6,4);
T{6} = strings(6,4);
T{7} = strings(6,4);

for i=1:4 % loop through AR,MR,AS,MS
    pathStr = NAMES{i};
    for j=0:5
        if i>2 && j==4
            continue
        end
        
        if j==5
            T{1}(j+1,i) = stats.N.(pathStr).GEQ.grade1.any;
            T{2}(j+1,i) = stats.N.(pathStr).GEQ.grade1.audMur;
            T{3}(j+1,i) = stats.N.(pathStr).GEQ.grade1.noAudMur;
            T{4}(j+1,i) = round(stats.P.audMur.(pathStr).GEQ.grade1,2);
            T{5}(j+1,i) = round(stats.P.audMur.(pathStr).GEQ.grade1,2);
            T{6}(j+1,i) = sprintf('%g/%g',T{2}(j+1,i),T{1}(j+1,i));
            T{7}(j+1,i) = sprintf('%g/%g',T{3}(j+1,i),T{1}(j+1,i));
            continue
        end
        gradeStr = strcat('grade',sprintf('%g',j));
        T{1}(j+1,i) = stats.N.(pathStr).EQ.(gradeStr).any;
        T{2}(j+1,i) = stats.N.(pathStr).EQ.(gradeStr).audMur;
        T{3}(j+1,i) = stats.N.(pathStr).EQ.(gradeStr).noAudMur;
        T{4}(j+1,i) = round(stats.P.audMur.(pathStr).EQ.(gradeStr),2);
        T{5}(j+1,i) = round(stats.P.(pathStr).EQ.(gradeStr).audMur,2);
        T{6}(j+1,i) = sprintf('%g/%g',T{2}(j+1,i),T{1}(j+1,i));
        T{7}(j+1,i) = sprintf('%g/%g',T{3}(j+1,i),T{1}(j+1,i));
        
    end
end
RnamesN = {'g0','g1','g2','g3','g4','presence'};
RnamesP = {'g0','g1','g2','g3','g4','presence'};
Cnames = {'AR','MR','AS','MS'};
T{1} = array2table(T{1},'VariableNames',Cnames,'RowNames',RnamesN);
T{2} = array2table(T{2},'VariableNames',Cnames,'RowNames',RnamesN);
T{3} = array2table(T{3},'VariableNames',Cnames,'RowNames',RnamesN);
T{4} = array2table(T{4},'VariableNames',Cnames,'RowNames',RnamesP);
T{5} = array2table(T{5},'VariableNames',Cnames,'RowNames',RnamesP);
T{6} = array2table(T{6},'VariableNames',Cnames,'RowNames',RnamesP);
T{7} = array2table(T{7},'VariableNames',Cnames,'RowNames',RnamesP);

T{1} = giveTitle2table(T{1},{'# VHD cases for different grades'});
T{2} = giveTitle2table(T{2},{'Number of VHD cases where murmur is present'});
T{3} = giveTitle2table(T{3},{'Number of VHD cases where murmur is NOT present'});
T{4} = giveTitle2table(T{4},{'fraction of VHD cases where murmur is present'});
T{5} = giveTitle2table(T{5},{'fraction of audible murmurs where VHD is present'});
T{6} = giveTitle2table(T{6},{'Number of murmurs relative to Number of VHD cases'});
T{7} = giveTitle2table(T{7},{'Number of inaudible murmurs relative to Number of VHD cases'});


stats.table.(genvarname('Number of VHD cases for different grades'))         = T{1};
stats.table.(genvarname('Number of VHD cases where murmur is present'))      = T{2};
stats.table.(genvarname('Number of VHD cases where murmur is NOT present'))  = T{3};
stats.table.(genvarname('fraction of VHD cases where murmur is present'))    = T{4};
stats.table.(genvarname('fraction of weakly audible murmurs where VHD is present')) = T{5};
stats.table.(genvarname('Number of murmurs relative to Number of VHD cases')) = T{6};
stats.table.(genvarname('Number of inaudible murmurs relative to Number of VHD cases')) = T{7};

stats.table.VHDprevalenceAndProbMur = [stats.table.NumberOfVHDCasesForDifferentGrades,...
stats.table.fractionOfVHDCasesWhereMurmurIsPresent] %#ok<*NOPTS>
stats.table.VHDprevalenceVHDAndVHDbutNoMur = [stats.table.NumberOfVHDCasesForDifferentGrades,...
stats.table.NumberOfVHDCasesWhereMurmurIsNOTPresent]
%% FOR WHICH VHD CAN DIASTOLIC MURMURS BE HEARD?
NAMES = {'AR','MR','AS','MS','total'};
T = strings(2,5);


T(1,1) = sum(and(stats.I.AR.GEQ.grade3.any,stats.I.diaMur.agree.any));
T(2,1) = sum(and(stats.I.AR.GEQ.grade3.any,stats.I.diaMur.obs1or2.any));
T(1,2) = sum(and(stats.I.MR.GEQ.grade3.any,stats.I.diaMur.agree.any));
T(2,2) = sum(and(stats.I.MR.GEQ.grade3.any,stats.I.diaMur.obs1or2.any));
T(1,3) = sum(and(stats.I.AS.GEQ.grade3.any,stats.I.diaMur.agree.any));
T(2,3) = sum(and(stats.I.AS.GEQ.grade3.any,stats.I.diaMur.obs1or2.any));
T(1,4) = sum(and(stats.I.MS.GEQ.grade3.any,stats.I.diaMur.agree.any));
T(2,4) = sum(and(stats.I.MS.GEQ.grade3.any,stats.I.diaMur.obs1or2.any));
T(1,5) = stats.N.diamur.agree.any;
T(2,5) = stats.N.diamur.obs1or2.any;
stats.table.diastolicSummary = array2table(T,'VariableNames',NAMES,'RowNames',{'agreed diast. murmur','obs. 1 or 2 heard diastolic murmur'})

%% GENERATE TABLE THAT SHOWS OVERLAP OF SIGNIFICANT VHS
% DEFINITIONS
% AR -- clinically interesting if GEQ 3
% MR -- clinically interesting if GEQ 3
% AS -- clinically interesting if GEQ 1 or higher (any degree is of interest)
% MS -- clinically interesting if GEQ 1 or higher (any degree is of interest)
NAMES = {'AR','MR','AS','MS'};
T = strings(4,4);

a = sum(and(stats.I.AR.GEQ.grade3.any,stats.I.AS.GEQ.grade1.any))
b = sum(or(stats.I.AR.GEQ.grade3.any,stats.I.AS.GEQ.grade1.any))

for i=1:4
    pathStr1 = NAMES{i};
    if pathStr1(2)=='R'
        gradeStr1 = 'grade3';
    else
        gradeStr1 = 'grade1';
    end
    
    for j=1:4
        if j<i
            T(i,j) = '--'
            continue
        end
        pathStr2 = NAMES{j};
        if pathStr2(2)=='R'
            gradeStr2 = 'grade3';
        else
            gradeStr2 = 'grade1';
        end
        
        pathStr2 = NAMES{j};
        a = sum(and(stats.I.(pathStr1).GEQ.(gradeStr1).any, ...
                    stats.I.(pathStr2).GEQ.(gradeStr2).any))
        b = sum(or(stats.I.(pathStr1).GEQ.(gradeStr1).any, ...
                   stats.I.(pathStr2).GEQ.(gradeStr2).any))
        str = sprintf('%g/%g',a,b);
        T(i,j) = str;
    end
end

Cnames = {'AR','MR','AS','MS'}; 
Rnames = Cnames;
stats.table.pathOverlap = array2table(T,'VariableNames',Cnames,'RowNames',Rnames);
stats.table.pathOverlap = giveTitle2table(stats.table.pathOverlap,{'overlap of clinically significant VHD'})
%% HOW IS PREDICTABILITY OF AR RELATED TO AUSCULTATION POSITION?
% Is murmur in one location more predictive of AR in some locations than in
% others?
clear T
VHDnames = {'AR','MR','AS','MS'};

for murThr=1:3
for k=1:4  
    vhd = VHDnames{k};
    
    numGrades = 4 - (k>2);
    
    T{1} = strings(numGrades,4);
    T{2} = nan(numGrades,2);
    T{3} = strings(numGrades,4);
    T{4} = strings(numGrades,4);
    T{5} = strings(numGrades,4);

    for aa=1:4
        for g=1:numGrades
        I1     = stats.I.(vhd).GEQ.(sprintf('grade%g',g)).any;
        murStr = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
        I2     = data.(murStr)>=murThr;
        
        [p,ci] = condProb(I1,I2,true);
        T{1}(g,aa) = displayEstWithCi(p, ci,"percentage","full", 2); 
        % compute murmur prevalence:
        T{2}(g,aa) = round(100*prob(I1), 1);
        % compute risk-ratio:
        [rr,ciRR] = relRisk(I1,I2,true);
        T{3}(g,aa) = displayEstWithCi(rr, ciRR,"decimal","full", 2)  
        
        % sensitivities and specificities:
        [p,ci] = condProb(I2,I1,true);
        T{4}(g,aa) = displayEstWithCi(p, ci,"percentage","full", 1);
        [p,ci] = condProb(~I2,~I1,true);
        T{5}(g,aa) = displayEstWithCi(p, ci,"percentage","full", 1);     
        end
    end
    
    % get number of cases in general population:
    for i=1:numGrades
        T{2}(i,2) = stats.N.(vhd).GEQ.(sprintf('grade%g',i)).any;
    end
    T{2} = T{2}(:,1:2);
    
    Rnames = {'grade 1','grade 2','grade 3','grade 4'}; 
    Cnames = {'loc 1','loc 2','loc 3','loc 4'};
    if numGrades==3
        Rnames = Rnames(1:3);
    end
    
    T{1} = array2table(T{1},'VariableNames',Cnames,'RowNames',Rnames);
    T{2} = array2table(T{2},'VariableNames',{'prob.','N'},'RowNames',Rnames);
    T{3} = array2table(T{3},'VariableNames',Cnames,'RowNames',Rnames);
    T{4} = array2table(T{4},'VariableNames',Cnames,'RowNames',Rnames);
    T{5} = array2table(T{5},'VariableNames',Cnames,'RowNames',Rnames);
    
    T{1} =  giveTitle2table(T{1},{sprintf('prob(%s grade >=i | max-murmur-grade >= %g)',vhd,murThr)});
    T{2} =  giveTitle2table(T{2},{sprintf('prevalence %s grade >=i',vhd)});
    T{3} =  giveTitle2table(T{3},{sprintf('Risk ratio of %s wrt. max-murmur-grade',vhd)});
    T{4} =  giveTitle2table(T{4},{sprintf('prob(max-murmur-grade >= %g| %s grade >=i )',murThr,vhd)});
    T{5} =  giveTitle2table(T{5},{sprintf('prob(max-murmur-grade < %g| %s grade <i )',murThr,vhd)});

    stats.table.(vhd).murmurRelationship.allAA.(num2name(murThr)) = [T{1},T{2},T{3},T{4},T{5}];
end
end

%% HOW IS PREDICTABILITY OF VHD RELATED TO maximum murmur grade?
% Is maximum murmur grade (maximum over all 4 auscultation locations)
% predictive of VHD?
clear T
VHDnames = {'AR','MR','AS','MS'};

for murThr=1:3
for k=1:4 % runs over VHD:s
    vhd = VHDnames{k};
    numGrades = 4 - (k>2);
    
    T{1} = strings(numGrades,1);
    T{2} = nan(numGrades,2);
    T{3} = strings(numGrades,1);
    T{4} = strings(numGrades,1);
    T{5} = strings(numGrades,1);

    for g=1:numGrades % runs over disease grades
        % I1 is disease, I2 is disease detection
        I1     = stats.I.(vhd).GEQ.(sprintf('grade%g',g)).any;
        I2     = data.maxMeanMurGrade>=murThr;
        
        [p,ci] = condProb(I1,I2,true);
        T{1}(g) = displayEstWithCi(p, ci,"percentage","full", 2); 
        % compute murmur prevalence:
        T{2}(g,1) = round(100*prob(I1), 1);
        % compute risk-ratio:
        [rr,ciRR] = relRisk(I1,I2,true);
        T{3}(g) = displayEstWithCi(rr, ciRR,"decimal","full", 2)  
        
        % sensitivities and specificities:
        [p,ci] = condProb(I2,I1,true);
        T{4}(g) = displayEstWithCi(p, ci,"percentage","full", 1);
        [p,ci] = condProb(~I2,~I1,true);
        T{5}(g) = displayEstWithCi(p, ci,"percentage","full", 1);
    end
    
    % get number of cases in general population:
    for i=1:numGrades
        T{2}(i,2) = stats.N.(vhd).GEQ.(sprintf('grade%g',i)).any;
    end
    
    T{2} = T{2}(:,1:2);
    
    Rnames = {'grade 1','grade 2','grade 3','grade 4'}; 
    Cnames = {' '};
    if numGrades==3
        Rnames = Rnames(1:3);
    end
    
    T{1} = array2table(T{1},'VariableNames',Cnames,'RowNames',Rnames);
    T{2} = array2table(T{2},'VariableNames',{'prob.','N'},'RowNames',Rnames);
    T{3} = array2table(T{3},'VariableNames',Cnames,'RowNames',Rnames);
    T{4} = array2table(T{4},'VariableNames',Cnames,'RowNames',Rnames);
    T{5} = array2table(T{5},'VariableNames',Cnames,'RowNames',Rnames);
    
    T{1} =  giveTitle2table(T{1},{sprintf('prob(%s grade >=i | max-murmur-grade >= %g)',vhd,murThr)});
    T{2} =  giveTitle2table(T{2},{sprintf('prevalence %s grade >=i',vhd)});
    T{3} =  giveTitle2table(T{3},{sprintf('Risk ratio of %s wrt. max-murmur-grade',vhd)});
    T{4} =  giveTitle2table(T{4},{sprintf('prob(max-murmur-grade >= %g| %s grade >=i )',murThr,vhd)});
    T{5} =  giveTitle2table(T{5},{sprintf('prob(max-murmur-grade < %g| %s grade <i )',murThr,vhd)});

    stats.table.(vhd).murmurRelationship.maxOnly.(num2name(murThr)) = [T{1},T{2},T{3},T{4},T{5}];
end
end