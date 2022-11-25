load('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');
HMMpar.Bmatrix = Springer_B_matrix;
HMMpar.piVector = Springer_pi_vector;
HMMpar.totObsDist = Springer_total_obs_distribution;
%% Get PCG Features:
data = HSdataTrain;


fs = 44100;
% define which set to explore;
set = data(stats.I.AS.GEQ.grade1.any,:);
k = k+1;
id = set.UNIKT_LOPENR(k);
% index AuscultationArea
aa = 1; 
x = wav2TS(id,aa);


[assignedStates,heartPar] = runSpringerSegmentationAlgorithm(x, fs, HMMpar,[]);

% get the features used in segmentation algorithm
% % [PCG_Features, featuresFs] = getSpringerPCGFeatures(x, fs);

% get distribution parameters; heart-rate and systole duration:
[heartPar.rate, heartPar.sysDuration] = getHeartRateSchmidt(x, fs);

[delta, ~, qt] = viterbiDecodePCG_Springer(PCG_Features, ...
    Springer_pi_vector, Springer_B_matrix, ...
    Springer_total_obs_distribution, heartRate, systolicTimeInterval, featuresFs);

max(delta(end,:))
assignedStates = expand_qt(qt, featuresFs, fs, length(x));
% var(abs(audioData(assignedStates==2)))/var(abs(audioData(assignedStates==1)))


% if(true)
%    figure('Name','Derived state sequence');
%    t1 = (1:length(audioData))./Fs;
%    plot(t1,normalise_signal(audioData),'k');
%    hold on;
%    plot(t1,assigned_states,'r--');
%    legend('Audio data', 'Derived states');
% end

% plot systole and diastole
clf
Nds = 20;
audioDataDS = downsample(x,Nds);
assignedStatesDS = downsample(assignedStates,Nds);
% subplot(211)
%     M = getScaleogram(audioDataDS,Fs,true);
%     hold on
%     Jsys = find(assigned_states1==2);
%     Jdia = find(assigned_states1==4);
%     scatter(Jsys/Nds,ones(1,length(Jsys))*40,'r.')
%     scatter(Jdia/Nds,ones(1,length(Jdia))*40,'g.')
% plot systole and diastole
% subplot(212)

audioDataDS = downsample(x,Nds);
assignedStatesDS = downsample(assignedStates,Nds);
mean(abs(audioDataDS(assignedStatesDS==1)))

% M = getScaleogram(audioDataDS,Fs,true);
% M = M(1:60,:);
% [cycleIndex,segLines] = states2cycles(assignedStatesDS);
% N.ds.all = 4; N.ds.initApprox = 10; N.smoothe = 2; N.rand = 50;
% xx{1} = segLines(1):segLines(2);
% xx{2} = segLines(2):segLines(3);
% Rsqr = maxModMatchScoreAll2d(segLines, M,1,N,false);
% round(mean(Rsqr),2)*100
close all
figure()
    getScaleogram(audioDataDS,fs,true);
    hold on
    Jsys = find(assignedStates==2);
    Jdia = find(assignedStates==4);
    scatter(Jsys/Nds,ones(1,length(Jsys))*40,'r.')
    scatter(Jdia/Nds,ones(1,length(Jdia))*40,'g.')