% In this script, I estimate the optimal decision threshold based on each
% training set that was used during cross validation.
%% load network, activiation thresholds, and segmentation info:
load networksCVnoNoiseDropout
load performanceCVnoNoiseDropOut
load nodes
%% show contents:
whos -file networksCVnoNoiseDropout
whos -file performanceCVnoNoiseDropOut
whos -file nodes
%% get training set performance and 
data0 = HSdataTrainAndVal;
Jnonoise = data0.MURMUR_1NOISE_REF_T72==0;
data0 = data0(Jnonoise,:);
nodes0.loc      = nodes.loc(Jnonoise,:);
nodes0.state    = nodes.state(Jnonoise,:);
nodes0.seglines = nodes.seglines(Jnonoise,:);

% define the training labels
aa = 1;
mg = 2;
pathThr = 2;
murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = data0.(murSt)>=mg;
N0 = height(data0);
Nsplits = 8;
Nval = floor(N0/Nsplits);
Jall = 1:N0;

nn.AS = zeros(1,Nsplits);
nn.AR = zeros(1,Nsplits);
nn.MS = zeros(1,Nsplits);
nn.MR = zeros(1,Nsplits);

alg.sensVec.AR = zeros(1,Nsplits);
alg.specVec.AR = zeros(1,Nsplits);
ann.sensVec.AR = zeros(1,Nsplits);
ann.specVec.AR = zeros(1,Nsplits);

for i=1:Nsplits
    i %#ok<*NOPTS>
    Jval   = (i-1)*Nval+1:i*Nval;
    dataVal   = data0(Jval,:);
    segLinesVal   = nodes0.seglines(Jval,:);

    N.ds = 20;
    N.cs = 4;
    N.os = 2;
    N.segPerPCG = 6;

    IposCasesVal = IposCases(Jval);
    Ncomp = [30,200];
    % get training set and validation set
    [Xval,Yval] = genTrainOrValSet(dataVal,N,aa,IposCasesVal,balanceVal,...
                                segLinesVal,Ncomp);
                            
    % get activations:
    Nloc = 1;
    showPerfCurve = true;
    Yval = data0.ARGRADE_T72(Jval)>=pathThr;
    [~,~,~,~,p] = performanceSummaryNeurNet(Xval,Yval,networks{i},...
                            Nval,N.segPerPCG,showPerfCurve,Nloc);
    % convert to logical:
    Yval = data0.ARGRADE_T72(Jval)>=pathThr;
    % get predictions using above estimate of optimal threshold:
    Ypred = p.murWhole>=uVec(i)*1.1;
    % estimate sensitivity and specificity:
    alg.sensVec.AR(i) = condProb(Ypred,Yval)
    alg.specVec.AR(i) = condProb(~Ypred,~Yval)
    ann.sensVec.AR(i) = condProb(data0.Murmur_1_grade_ref_ny_T72(Jval)>=mg,Yval)
    ann.specVec.AR(i) = condProb(data0.Murmur_1_grade_ref_ny_T72(Jval)<mg,~Yval)
end