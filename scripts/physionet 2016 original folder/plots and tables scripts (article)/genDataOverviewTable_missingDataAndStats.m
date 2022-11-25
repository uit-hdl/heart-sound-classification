% in this script I generate the data used in the table in the appendix that
% provides an overwiev of the dataset.
%% PREVALENCE OF VHD 
data = HSdata(~HSdata.noiseOnly,:);
% how many with AR for each grade?
x = []
for i=1:4
    x = [x,data.MSgrade>=i];
end
sum(x)
100*round(mean(x),4)
%%
diseaseNames = ["AR" "MR" "AS" "MS"];
VHDthr = 2;
for i=1:4
    VHDname = diseaseNames(i);
    varName = sprintf('%sgrade',VHDname);
    sprintf('%g (%g%%)',round(sum(HSdata.(varName)==VHDthr),2),...
                    dec2perc(mean(HSdata.(varName)==VHDthr),2))
end
%% what percentage of AR were female?
[p,ci] = condProb(HSdata.sex==1,HSdata.MSgrade>=1,true)
dec2perc(p,2)
round(ci2SD(ci)*100,1)
%% How many significant cases were male and female?
varName = 'MSgrade';
sex = 0;
thr = (varName(2)=="R")*3 + (varName(2)=="S")*1;
sum(and(HSdata.sex==sex,HSdata.(varName)>=thr))
dec2perc(condProb(HSdata.sex==sex,HSdata.(varName)>=1,true),1)
%% age distribution
computeCImeanEst(HSdata.age(HSdata.MSgrade>=1),"3")

% quantile(HSdata.age(HSdata.ARgrade>=3),[1/3,2/3])
IageG{1} = and(HSdata.age>=40,HSdata.age<50)
IageG{2} = and(HSdata.age>=50,HSdata.age<60)
IageG{3} = and(HSdata.age>=60,HSdata.age<70)
IageG{4} = and(HSdata.age>=70,HSdata.age<80)
IageG{5} = HSdata.age>=80

ageG = IageG{1};
VHDname = 'MR';
varName = sprintf('%sgrade',VHDname);
VHDthr = 3;

for i=1:5
    sprintf('%g (%g%%)',sum(and(HSdata.(varName)>=VHDthr,IageG{i})),...
                    dec2perc(condProb(IageG{i},HSdata.(varName)>=VHDthr),1))
end
%% presence of symptoms, all participants
X.dyspneaCalmlyFlat = dec2perc(mean(isnan(HSdata.dyspneaCalmlyFlat)),2)
X.dyspneaRest = dec2perc(mean(isnan(HSdata.dyspneaRest)),2)
X.dyspneaRestOrFlat = dec2perc(mean(isnan(HSdata.dyspneaRestOrFlat)),2)
X.highBP = dec2perc(mean(isnan(HSdata.highBP)),2)
X.angina = dec2perc(mean(isnan(HSdata.angina)),2)
X.pulseSpiro = dec2perc(mean(isnan(HSdata.pulseSpiro)),2)
X.diabetes = dec2perc(mean(isnan(HSdata.diabetes)),2)
X.highBP = dec2perc(mean(isnan(HSdata.highBP)),2)
X.chestPain = dec2perc(mean(isnan(HSdata.chestPain)),2)
X.smoke = dec2perc(mean(isnan(HSdata.smoke)),2)

T = struct2table(X)


%% presence of symptoms, VHD-subgroups
clear X

clinicVar = "smoke";
clinicVar = "dyspneaRestOrFlat";
diseaseNames = ["AR" "MR" "AS" "MS"];
data = HSdata(~isnan(HSdata.(clinicVar)),:);

for i=1:4
    var = diseaseNames(i);
    varName = sprintf('%sgrade',var);
    thr = (varName(2)=="R")*3 + (varName(2)=="S")*1;

    % n confirmed dyspnea in positive class
    X.(var).NconfirmedSymp = sum(data.(clinicVar)(data.(varName)>=thr));
    X.(var).PropConfirmedSymp = dec2perc(mean(data.(clinicVar)(data.(varName)>=thr)),1);
    % n missing data in positive class
    X.(var).Nmissing = sum(and(isnan(HSdata.(clinicVar)),HSdata.(varName)>=thr))
    % proportion missing in positive class
    X.(var).propMissing = dec2perc(condProb(...
        isnan(HSdata.(clinicVar)),HSdata.(varName)>=thr) ,1)
    X.(var).presentMissing = sprintf('%g (%g%%)',X.(var).Nmissing,X.(var).propMissing)
    X.(var).presentSympt = sprintf('%g (%g%%)',X.(var).NconfirmedSymp,X.(var).PropConfirmedSymp)
end

sprintf('%s %s %s %s',X.AR.presentSympt,X.MR.presentSympt,X.AS.presentSympt,X.MS.presentSympt)
sprintf('%s %s %s %s',X.AR.presentMissing,X.MR.presentMissing,X.AS.presentMissing,X.MS.presentMissing)

%% presence of symptoms, continuous variables
clear X

clinicVar = "costumVar";
HSdata.costumVar = HSdata.murGrade4>0;
% clinicVar = "dyspneaRestOrFlat";
diseaseNames = ["AR" "MR" "AS" "MS"];
data = HSdata(~isnan(HSdata.(clinicVar)),:);

for i=1:4
    var = diseaseNames(i);
    varName = sprintf('%sgrade',var);
    thr = (varName(2)=="R")*3 + (varName(2)=="S")*1;

    % mean
    X.(var).mean = round(mean(data.(clinicVar)(data.(varName)>=thr)),1);
    ci = computeCImeanEst(data.(clinicVar)(data.(varName)>=thr));
    X.(var).meanSD = round(ci2SD(ci),1);
    
    X.(var).Nmissing = sum(and(isnan(HSdata.(clinicVar)),HSdata.(varName)>=thr))
    % proportion missing in positive class
    X.(var).propMissing = dec2perc(condProb(...
        isnan(HSdata.(clinicVar)),HSdata.(varName)>=thr) ,1)
    
    X.(var).presentMissing = sprintf('%g (%g%%)',X.(var).Nmissing, X.(var).propMissing)
    X.(var).presentSympt = sprintf('%g (%g)',X.(var).mean, X.(var).meanSD)
end

sprintf('%s %s %s %s',X.AR.presentSympt,X.MR.presentSympt,X.AS.presentSympt,X.MS.presentSympt)
sprintf('%s %s %s %s',X.AR.presentMissing,X.MR.presentMissing,X.AS.presentMissing,X.MS.presentMissing)



%% MURMUR PREVALENCE IN VHD's
murThr = 1;
vhdThr = 4;
% n murmurs for each VHD:
n = sum([and(HSdata.ARgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr),...
     and(HSdata.MRgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr),...
     and(HSdata.ASgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr),...
     and(HSdata.MSgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr)])'
 
 % prevalence murmur for each VHD:
p = dec2perc([condProb(HSdata.maxMeanMurGrade>=murThr,HSdata.ARgrade>=vhdThr),...
     condProb(HSdata.maxMeanMurGrade>=murThr,HSdata.MRgrade>=vhdThr),...
     condProb(HSdata.maxMeanMurGrade>=murThr,HSdata.ASgrade>=vhdThr),...
     condProb(HSdata.maxMeanMurGrade>=murThr,HSdata.MSgrade>=vhdThr)],1)'
 
 [n,p]
 
 sprintf('%g (%g%%)  %g (%g%%)  %g (%g%%)  %g (%g%%)',...
     n(1),p(1),n(2),p(2),n(3),p(3),n(4),p(4))

 
 % relative risk:
estCi = true;
vhdThr = 2;
[RR,RRci] = relRisk(HSdata.ARgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr,estCi);
RRci = round(RRci,2);
RR   = round(RR,2);
sprintf('%g (%g-%g)',RR,RRci(1),RRci(2))

[RR,RRci] = relRisk(HSdata.MRgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr,estCi);
RRci = round(RRci,2);
RR   = round(RR,2);
sprintf('%g (%g-%g)',RR,RRci(1),RRci(2))

[RR,RRci] = relRisk(HSdata.ASgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr,estCi);
RRci = round(RRci,2);
RR   = round(RR,2);
sprintf('%g (%g-%g)',RR,RRci(1),RRci(2))

[RR,RRci] = relRisk(HSdata.MSgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr,estCi);
RRci = round(RRci,2);
RR   = round(RR,2);
sprintf('%g (%g-%g)',RR,RRci(1),RRci(2))
 
 
 condProb(HSdata.MSgrade>=vhdThr,HSdata.maxMeanMurGrade>=murThr)/...
 condProb(HSdata.MSgrade>=vhdThr,~HSdata.maxMeanMurGrade>=murThr)

 
 
