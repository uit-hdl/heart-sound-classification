
load('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');
HMMpar.Bmatrix = Springer_B_matrix;
HMMpar.piVector = Springer_pi_vector;
HMMpar.totObsDist = Springer_total_obs_distribution;
data = HSdataTrain;
newk = 1;
k = 0;
colMap = parula(255);
%%
NdsAcf = 35;
Nds0 = 20;
Fs = 44100;
fs = floor(Fs/Nds0);
% define which set to explore;
set = data(stats.I.AS.GEQ.grade0.any,:);
redo = 0;
if newk==1
    k = k + 1 - redo;
end
id = set.UNIKT_LOPENR(k);
Nds  = floor(30/Nds0);

[states, fitInfo] = ngbrSegment(id,Nds0,NdsAcf,HMMpar,segNgbr,segNoNgbr,plt,Fs);


