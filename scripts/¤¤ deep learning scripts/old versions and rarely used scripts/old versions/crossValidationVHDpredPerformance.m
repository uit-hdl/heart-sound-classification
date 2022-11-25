%% initial
load nodes
load networks
clearvars Jval

data0 = HSdataTrainAndVal;
Jnonoise = data0.MURMUR_1NOISE_REF_T72==0;
data0 = data0(Jnonoise,:);
nodes0.loc      = nodes.loc(Jnonoise,:);
nodes0.state    = nodes.state(Jnonoise,:);
nodes0.seglines = nodes.seglines(Jnonoise,:);
%%
% define the labels
aa = 1;
mg = 2;
murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = data0.(murSt)>=mg;
N0 = height(data0);
NvalSets = 9;
Nval = floor(N0/NvalSets);
Jall = 1:N0;

Xval = cell(NvalSets);
Yval = cell(NvalSets);

% set to true if you want training sets instead of validation sets:
getTraining = false;
for i=1:NvalSets
    i %#ok<NOPTS>
    %%% get validation set %%%
    Jval      = (i-1)*Nval+1:i*Nval;
    if getTraining
        Jval      = setdiff(Jall,Jval);
    end
    dataVal   = data0(Jval,:);
    segLinesVal = nodes0.seglines(Jval,:);

    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;

    IposCasesVal = IposCases(Jval);
    balanceVal   = false;
    Ncomp = [30,200];

    [Xval{i},Yval{i}] = genTrainOrValSet(dataVal,N,aa,IposCasesVal,balanceVal,...
                                segLinesVal,Ncomp);
end

%%  check that everything is in order:
u_decision = zeros(1,NvalSets);
NpcgPerLoc = numel(Jval);
for i=1:NvalSets
    i
Nloc = 1;
figure()
[AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval{i},Yval{i},...
                          networks{i},NpcgPerLoc,N.segPerPCG,true,Nloc);
title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
% predict optimal decision threshold:

% get index where optimal tradeoff is attained (given sensitivity >= .85):
minSensitivity = 0.85;
[~,Ioptimal] = max( (1-X.whole + Y.whole).*(Y.whole>=minSensitivity));
% get decision optimal threshold:
u_decision(i)     = T.whole(Ioptimal) %#ok<NOPTS>
Ypred = median(p.murWhole>u_decision,2);

Nloc = 1;
Yactual  = mode(reshapeSegmentScores(Yval{i},NpcgPerLoc,N.segPerPCG,Nloc),2);

[condProb(Ypred,Yactual),condProb(~Ypred,~Yactual)]
end

%% predict diseases:

% storage cell for activations of each validation set:
activations = cell(1,NvalSets);
NpcgPerLoc = numel(Jval);
for i=1:NvalSets
figure()
[AUC,X,Y,T,p] = performanceSummaryNeurNet(Xval{i},Yval{i},...
                            networks{i},NpcgPerLoc,N.segPerPCG,true,Nloc);
title(sprintf('AUCseg=%.3g, AUCwhole=%.3g, location=%g',AUC.seg,AUC.whole,aa))
activations{i} = p.murWhole;
end






