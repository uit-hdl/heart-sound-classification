function [logLik,HeartPar,assignedStates] = tryAllHeartParameters(id,springerPar,Fs)
% estimates heart rate and systole length for all 4 auscultation areas. For
% each set of heart-sound parameters and for each of the 4 recording,
% estimates the state vector. Also saves the logLikelihood for each fit.

% Fs = 44100;
% Ind = stats.I.AR.GEQ.grade4.any;
% set = HSdataTrain(Ind,:);
% k=k+1;
% id = set.UNIKT_LOPENR(k);

HeartPar.sysLen    = zeros(4,1);
HeartPar.heartRate = zeros(4,1);
HeartPar.LogLik    = zeros(4,1);
% PCG_Features = cell(4,1);
% featuresFs   = cell(4,1);
assignedStates = cell(4,4);
logLik          = zeros(4,4);
for i=1:4
    audioData = wav2TS(id,i);
    [heartRate, systolicTimeInterval] = getHeartRateSchmidt(audioData, Fs);
    HeartPar.heartRate(i) = heartRate;
    HeartPar.sysLen(i) = systolicTimeInterval;
    
    for j=1:4
        audioData = wav2TS(id,j);
        [PCG_Features, featuresFs] = getSpringerPCGFeatures(audioData, Fs);
        [delta, ~, qt] = viterbiDecodePCG_Springer(PCG_Features, springerPar.piVector,...
                                           springerPar.Bmatrix, springerPar.totObsDist,...
                                           heartRate, systolicTimeInterval,...
                                           featuresFs);
        logLik(i,j) = max(delta(end,:));                         
        assignedStates{i,j} = expand_qt(qt, featuresFs, Fs, length(audioData));
    end
end
end


% for i=1:4
%     audioData = wav2TS(id,i);
%     [~, ~, qt] = viterbiDecodePCG_Springer(PCG_Features{i}, springerPar.piVector,...
%                                            springerPar.Bmatrix, springerPar.totObsDist,...
%                                            HeartPar.heartRate(i), HeartPar.sysLen(i),...
%                                            featuresFs{i});
%                                        
%     assignedStates{i} = expand_qt(qt, featuresFs{i}, Fs, length(audioData));
%     
% end




