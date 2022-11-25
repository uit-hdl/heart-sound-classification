function [DELTA,bestAA,worstAA,heartRate,systolicTimeInterval] = ...
                                        findBestAuscArea(id,springerPar,Fs,compDelta)
%% description
% computes the log likelihood of the estimated state sequence of each
% recording, and uses this to estimate which auscultation position is best
% to use. 

% Note: not in use.
%% preliminary
if nargin==2
    Fs = 44100;
    compDelta = false;
elseif nargin==3
    compDelta = false;
end

DELTA = nan(4,1);
heartRate = zeros(4,1);
systolicTimeInterval = zeros(4,1);
for aa=1:4
    audioData = wav2TS(id,aa);

    [PCG_Features, featuresFs] = getSpringerPCGFeatures(audioData, Fs);
    
    % Get PCG heart rate
    [heartRate(aa), systolicTimeInterval(aa)] = getHeartRateSchmidt(audioData, Fs);
    if compDelta
        [delta, ~, ~] = viterbiDecodePCG_Springer(PCG_Features, springerPar.piVector, ...
                         springerPar.Bmatrix, springerPar.totObsDist, ...
                         heartRate(aa), systolicTimeInterval(aa), featuresFs);
        DELTA(aa) = max(delta(end,:));
    end
end
[~,bestAA] = max(DELTA);
[~,worstAA] = min(DELTA);
end