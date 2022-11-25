% when doing cross validation, the ninth training set turned out to
% consistently deliver very poor results; not validation results only, but
% training results, and in all likelyhood the poor validation results were
% due to to poorly trained algorithm, and the issue lies in the training
% set, not the validation set. I take a closer look here to figure out what
% is special about this data set. Does it have fewer positive cases? How do
% the positive cases look like? Are some of them particularly nasty? this
% is the most likley factor to cause issues during training, since I
% duplicate the positive samples in order to balance the training set.
data0 = HSdataTrainAndVal;
% define the labels
aa = 1;
mg = 2;
murSt = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
IposCases = data0.(murSt)>=mg;
N0 = height(data0);
NvalSets = 9;
Nval = floor(N0/NvalSets);
Jall = 1:N0;



i = 9;
Jval   = (i-1)*Nval+1:i*Nval;
Jtrain = setdiff(Jall,Jval);
dataTrain = data0(Jtrain,:);
dataVal   = data0(Jval,:);
segLinesTrain = nodes.seglines(Jtrain,:);
segLinesVal   = nodes.seglines(Jval,:);

[sum(dataTrain.Murmur_1_grade_ref_ny_T72>=2),...
sum(dataTrain.MURMUR_1NOISE_REF_T72>0),...
sum(and(dataTrain.MURMUR_1NOISE_REF_T72>0,...
    dataTrain.Murmur_1_grade_ref_ny_T72>=2)),...
    sum(dataTrain.ASGRADE_T72>0),...
    mean(dataTrain.Murmur_1_grade_ref_ny_T72(dataTrain.Murmur_1_grade_ref_ny_T72>=2)),...
    ]
% training set 9 summary:
% 105 murmurs, 101 noisy samples




