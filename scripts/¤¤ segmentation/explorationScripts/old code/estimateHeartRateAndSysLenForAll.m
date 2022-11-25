
load('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');
springerPar.Bmatrix = Springer_B_matrix;
springerPar.piVector = Springer_pi_vector;
springerPar.totObsDist = Springer_total_obs_distribution;

n = height(HSdataTrain);
Fs = 44100;
% stores he estimated systolic lengths for each subject. Each entry
% contains a vector of size 4 that stores values for each location.
HeartPar.sysLen    = zeros(n,4);
HeartPar.heartRate = zeros(n,4);
HeartPar.LogLik    = zeros(n,4);

for i=1:n
    id = HSdataTrain.UNIKT_LOPENR(i);
    [~,~,~,heartRate,systolicTimeInterval] = findBestAuscArea(id,springerPar,Fs,false);
    HeartPar.HeartRate(i,:) = heartRate;
    HeartPar.sysLen(i,:) = systolicTimeInterval;
    HeartPar.LogLik(i,:) = DELTA;
    i
end

