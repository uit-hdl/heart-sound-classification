% Create table providing overwiev over annotator-murmur grade AUC for
% predicting AS, and effect on prediction of taking mean and max.
%#ok<*NOPTS>
clear AUCstruct
classThr = 1;
Ypred = HSdata.maxMeanMurGrade;
Ytarget = HSdata.ASgrade>=classThr;

% *** SA predictions ***
SAmax = zeros(height(HSdata),4);
for aa=1:5
    if aa<5
        murStr = sprintf('MURMUR_%gGRADENRSA_T72',aa);
        Inoise = HSdata.(sprintf('noise%g',aa));
        Ypred = HSdata.(murStr)(~Inoise);
        Ytarget = HSdata.ASgrade(~Inoise)>=classThr;
        AUC = performanceSummaryNeurNet([],Ytarget,Ypred,...
                                numel(Ypred),[],plotTrain);
        AUCstruct.SA(aa,1) = round(AUC.whole,3);
        SAmax(~Inoise,aa) = Ypred;
    else
        Inoise = HSdata.noiseOnly;
        SAmax = max(SAmax,[],2);
        Ypred = SAmax(~Inoise);
        Ytarget = HSdata.ASgrade(~Inoise)>=classThr;
        [AUC] = performanceSummaryNeurNet([],Ytarget,Ypred,...
                             numel(Ypred),[],plotTrain);
        AUCstruct.SA(aa,1) = round(AUC.whole,2);
    end
end

% *** AD predictions ***
ADmax = zeros(height(HSdata),4);
for aa=1:5
    if aa<5
        murStr = sprintf('MURMUR_%gGRADENRAD_T72',aa);
        Inoise = HSdata.(sprintf('noise%g',aa));
        Ypred = HSdata.(murStr)(~Inoise);
        Ytarget = HSdata.ASgrade(~Inoise)>=classThr;
        AUC = performanceSummaryNeurNet([],Ytarget,Ypred,...
                                numel(Ypred),[],plotTrain);
        AUCstruct.AD(aa,1) = round(AUC.whole,3);
        ADmax(~Inoise,aa) = Ypred;
    else
        Inoise = HSdata.noiseOnly;
        ADmax = max(ADmax,[],2);
        Ypred = ADmax(~Inoise);
        Ytarget = HSdata.ASgrade(~Inoise)>=classThr;
        [AUC] = performanceSummaryNeurNet([],Ytarget,Ypred,...
                             numel(Ypred),[],plotTrain);
        AUCstruct.AD(aa,1) = round(AUC.whole,2);
    end
end

% *** mean predictions ***
meanMax = zeros(height(HSdata),4);
for aa=1:5
    if aa<5
        murStr = sprintf('murGrade%g',aa);
        Inoise = HSdata.(sprintf('noise%g',aa));
        Ypred = HSdata.(murStr)(~Inoise);
        Ytarget = HSdata.ASgrade(~Inoise)>=classThr;
        [AUC] = performanceSummaryNeurNet([],Ytarget,Ypred,...
                                numel(Ypred),[],plotTrain);
        AUCstruct.mean(aa,1) = round(AUC.whole,3);
        meanMax(~Inoise,aa) = Ypred;
    else
        Inoise = HSdata.noiseOnly;
        meanMax = max(meanMax,[],2);
        Ypred = meanMax(~Inoise);
        Ytarget = HSdata.ASgrade(~Inoise)>=classThr;
        [AUC] = performanceSummaryNeurNet([],Ytarget,Ypred,...
        numel(Ypred),[],plotTrain);
        AUCstruct.mean(aa,1) = round(AUC.whole,2);
    end
end

% *** predicted mean-murmur-grade ***
AUCstruct.pred(1,1) = 96.3;
AUCstruct.pred(2,1) = 94.7;
AUCstruct.pred(3,1) = 97.3;
AUCstruct.pred(4,1) = 93.2;
AUCstruct.pred(5,1) = 97.9;
T = table(struct2table(AUCstruct),'v',{'overwiev of AUC (AS>=1) for different variables'},...
    'r',{'pos1','pos2','pos3','pos4','max'})

meanVsMaxAUC = 100*(AUCstruct.mean-max([AUCstruct.AD,AUCstruct.SA],[],2));
T = addvars(T,meanVsMaxAUC)
