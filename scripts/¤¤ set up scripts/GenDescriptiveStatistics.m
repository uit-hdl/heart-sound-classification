data = HSdata;
% how many have detactable murmur?
statsAll.I.audMurWeak.any = data.maxMeanMurGrade>=1;
statsAll.N.any = height(data);
statsAll.N.audMurWeak.any = sum(statsAll.I.audMurWeak.any);
definitions.audMurWeak = sprintf('The maximum of the mean murmur grades is greater than or equal to 1.\nNote that this is a very weak definition, and requires only one observer to hear murmur in location.');

% ROWS WITH DIASTOLIC MURMURS
IdiaMurObs1 = 0;
IdiaMurObs2 = 0;
for i=1:4
    str1 = sprintf('MURMUR_%gDIAAD_T72',i)
    str2 = sprintf('MURMUR_%gDIASA_T72',i)
    IdiaMurObs1 = IdiaMurObs1 + data.(str1);
    IdiaMurObs2 = IdiaMurObs2 + data.(str2);
end
statsAll.genData.maxDiaMurGrade = max(IdiaMurObs1,IdiaMurObs2);
statsAll.I.diaMur.Obs1.any = IdiaMurObs1>0;
statsAll.I.diaMur.Obs2.any = IdiaMurObs2>0;
statsAll.I.diaMur.agree.any  = min(statsAll.I.diaMur.Obs1.any,statsAll.I.diaMur.Obs2.any);
statsAll.I.diaMur.obs1or2.any = max(statsAll.I.diaMur.Obs1.any,statsAll.I.diaMur.Obs2.any);
statsAll.N.diamur.agree.any = sum(statsAll.I.diaMur.agree.any);
statsAll.N.diamur.obs1or2.any = sum(statsAll.I.diaMur.obs1or2.any);

% AR,MR,AS,MS indicator vectors and numbers for different grades
VHDnames = {'AR','MR','AS','MS'};

for i=1:4
    VHDStr  = VHDnames{i};
    statsAll.I.(VHDnames{i}).presence.any = data.(strcat(VHDnames{i},'PRESENCE_T72'));
    statsAll.N.(VHDnames{i}).presence.any = sum(statsAll.I.(VHDnames{i}).presence.any,'omitnan');
    % cycle through grades
    for j=0:4
        if i>2 && j==4
            continue
        end

        gradeStr= sprintf('grade%g',j);
        statsAll.I.(VHDStr).GEQ.(gradeStr).any = (data.(strcat(VHDStr,'grade'))>=j);
        statsAll.I.(VHDStr).EQ.(gradeStr).any  = (data.(strcat(VHDStr,'grade'))==j);
        % how many with pathology have grade of so-and-so?
        statsAll.N.(VHDStr).GEQ.(gradeStr).any = sum(statsAll.I.(VHDStr).GEQ.(gradeStr).any);
        statsAll.N.(VHDStr).EQ.(gradeStr).any  = sum(statsAll.I.(VHDStr).EQ.(gradeStr).any);
        
        % how many with PATH of GRADE (relation to) x had AUDIBLE MURMUR?
        statsAll.N.(VHDStr).GEQ.(gradeStr).audMurWeak = sum(and(statsAll.I.(VHDStr).GEQ.(gradeStr).any,...
                                                             statsAll.I.audMurWeak.any));
        statsAll.N.(VHDStr).GEQ.(gradeStr).noAudMurWeak = sum(and(statsAll.I.(VHDStr).GEQ.(gradeStr).any,...
                                                              ~statsAll.I.audMurWeak.any));
                                                          
        
        statsAll.N.(VHDStr).EQ.(gradeStr).audMurWeak = sum(and(statsAll.I.(VHDStr).EQ.(gradeStr).any,...
                                                             statsAll.I.audMurWeak.any));
        statsAll.N.(VHDStr).EQ.(gradeStr).noAudMurWeak = sum(and(statsAll.I.(VHDStr).EQ.(gradeStr).any,...
                                                           ~statsAll.I.audMurWeak.any));
        
        % P(PATH of GRADE >= x | AUDIBLE MURMUR)?                                       
        statsAll.P.(VHDStr).GEQ.(gradeStr).audMurWeak = statsAll.N.(VHDStr).GEQ.(gradeStr).audMurWeak/...
                                                   statsAll.N.audMurWeak.any ;
        % P(AUDIBLE MURMUR | PATH of GRADE >= x)?                               
        statsAll.P.audMurWeak.(VHDStr).GEQ.(gradeStr) = statsAll.N.(VHDStr).GEQ.(gradeStr).audMurWeak/...
                                                   statsAll.N.(VHDStr).GEQ.(gradeStr).any ;
        
        % P(PATH of GRADE == x | AUDIBLE MURMUR)?                                       
        statsAll.P.(VHDStr).EQ.(gradeStr).audMurWeak = statsAll.N.(VHDStr).EQ.(gradeStr).audMurWeak/...
                                                   statsAll.N.audMurWeak.any ;
        % P(AUDIBLE MURMUR | PATH of GRADE == x)?                               
        statsAll.P.audMurWeak.(VHDStr).EQ.(gradeStr) = statsAll.N.(VHDStr).EQ.(gradeStr).audMurWeak/...
                                                   statsAll.N.(VHDStr).EQ.(gradeStr).any ; 
                                                         
    end
end

definitions.red = sprintf('Refers to the data set which only contains rows for which echo results are available.');
definitions.GEQ = sprintf('Greater than or Eual to.');
%% CREATE DESCRIPTIVE TABLES
clear T
NAMES = {'AR','MR','AS','MS'};
T{1} = nan(6,4);
T{2} = nan(6,4);
T{3} = nan(6,4);
T{4} = nan(6,4);
T{5} = nan(6,4);
T{6} = strings(6,4);
T{7} = strings(6,4);

for i=1:4
    pathStr = NAMES{i};
    for j=0:5
        if i>2 && j==4
            continue
        end
        
        if j==5
            T{1}(j+1,i) = statsAll.N.(pathStr).GEQ.grade1.any;
            T{2}(j+1,i) = statsAll.N.(pathStr).GEQ.grade1.audMurWeak;
            T{3}(j+1,i) = statsAll.N.(pathStr).GEQ.grade1.noAudMurWeak;
            T{4}(j+1,i) = round(statsAll.P.audMurWeak.(pathStr).GEQ.grade1,2);
            T{5}(j+1,i) = round(statsAll.P.audMurWeak.(pathStr).GEQ.grade1,2);
            T{6}(j+1,i) = sprintf('%g/%g',T{2}(j+1,i),T{1}(j+1,i));
            T{7}(j+1,i) = sprintf('%g/%g',T{3}(j+1,i),T{1}(j+1,i));
            continue
        end
        gradeStr = strcat('grade',sprintf('%g',j));
        T{1}(j+1,i) = statsAll.N.(pathStr).EQ.(gradeStr).any;
        T{2}(j+1,i) = statsAll.N.(pathStr).EQ.(gradeStr).audMurWeak;
        T{3}(j+1,i) = statsAll.N.(pathStr).EQ.(gradeStr).noAudMurWeak;
        T{4}(j+1,i) = round(statsAll.P.audMurWeak.(pathStr).EQ.(gradeStr),2);
        T{5}(j+1,i) = round(statsAll.P.(pathStr).EQ.(gradeStr).audMurWeak,2);
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

T{1} = giveTitle2table(T{1},{'Number of VHD cases for different grades'});
T{2} = giveTitle2table(T{2},{'Number of VHD cases where murmur is present'});
T{3} = giveTitle2table(T{3},{'Number of VHD cases where murmur is NOT present'});
T{4} = giveTitle2table(T{4},{'fraction of VHD cases where murmur is present'});
T{5} = giveTitle2table(T{5},{'fraction of weakly audible murmurs where VHD is present'});
T{6} = giveTitle2table(T{6},{'Number of murmurs relative to Number of VHD cases'});
T{7} = giveTitle2table(T{7},{'Number of inaudible murmurs relative to Number of VHD cases'});


statsAll.table.(genvarname('Number of VHD cases for different grades'))         = T{1};
statsAll.table.(genvarname('Number of VHD cases where murmur is present'))      = T{2};
statsAll.table.(genvarname('Number of VHD cases where murmur is NOT present'))  = T{3};
statsAll.table.(genvarname('fraction of VHD cases where murmur is present'))    = T{4};
statsAll.table.(genvarname('fraction of weakly audible murmurs where VHD is present')) = T{5};
statsAll.table.(genvarname('Number of murmurs relative to Number of VHD cases')) = T{6};
statsAll.table.(genvarname('Number of inaudible murmurs relative to Number of VHD cases')) = T{7};

statsAll.table.VHDprevalenceAndProbMur = [statsAll.table.NumberOfVHDCasesForDifferentGrades,...
statsAll.table.fractionOfVHDCasesWhereMurmurIsPresent] %#ok<*NOPTS>
statsAll.table.VHDprevalenceVHDAndVHDbutNoMur = [statsAll.table.NumberOfVHDCasesForDifferentGrades,...
statsAll.table.NumberOfVHDCasesWhereMurmurIsNOTPresent]
%% FOR WHICH VHD CAN DIASTOLIC MURMURS BE HEARD?
NAMES = {'AR','MR','AS','MS','total'};
T = strings(2,5);

T(1,1) = sum(and(statsAll.I.AR.GEQ.grade3.any,statsAll.I.diaMur.agree.any));
T(2,1) = sum(and(statsAll.I.AR.GEQ.grade3.any,statsAll.I.diaMur.obs1or2.any));
T(1,2) = sum(and(statsAll.I.MR.GEQ.grade3.any,statsAll.I.diaMur.agree.any));
T(2,2) = sum(and(statsAll.I.MR.GEQ.grade3.any,statsAll.I.diaMur.obs1or2.any));
T(1,3) = sum(and(statsAll.I.AS.GEQ.grade3.any,statsAll.I.diaMur.agree.any));
T(2,3) = sum(and(statsAll.I.AS.GEQ.grade3.any,statsAll.I.diaMur.obs1or2.any));
T(1,4) = sum(and(statsAll.I.MS.GEQ.grade3.any,statsAll.I.diaMur.agree.any));
T(2,4) = sum(and(statsAll.I.MS.GEQ.grade3.any,statsAll.I.diaMur.obs1or2.any));
T(1,5) = statsAll.N.diamur.agree.any;
T(2,5) = statsAll.N.diamur.obs1or2.any;
statsAll.table.diastolicSummary = array2table(T,'VariableNames',NAMES,'RowNames',{'agreed diast. murmur','obs. 1 or 2 heard diastolic murmur'})

%% GENERATE TABLE THAT SHOWS OVERLAP OF SIGNIFICANT VHS
% DEFINITIONS
% AR -- clinically interesting if GEQ 3
% MR -- clinically interesting if GEQ 3
% AS -- clinically interesting if GEQ 1 or higher (any degree is of interest)
% MS -- clinically interesting if GEQ 1 or higher (any degree is of interest)
NAMES = {'AR','MR','AS','MS'};
T = strings(4,4);

a = sum(and(statsAll.I.AR.GEQ.grade3.any,statsAll.I.AS.GEQ.grade1.any))
b = sum(or(statsAll.I.AR.GEQ.grade3.any,statsAll.I.AS.GEQ.grade1.any))

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
        a = sum(and(statsAll.I.(pathStr1).GEQ.(gradeStr1).any, ...
                    statsAll.I.(pathStr2).GEQ.(gradeStr2).any))
        b = sum(or(statsAll.I.(pathStr1).GEQ.(gradeStr1).any, ...
                   statsAll.I.(pathStr2).GEQ.(gradeStr2).any))
        str = sprintf('%g/%g',a,b);
        T(i,j) = str;
    end
end

Cnames = {'AR','MR','AS','MS'}; 
Rnames = Cnames;
statsAll.table.pathOverlap = array2table(T,'VariableNames',Cnames,'RowNames',Rnames);
statsAll.table.pathOverlap = giveTitle2table(statsAll.table.pathOverlap,{'overlap of clinically significant VHD'})


