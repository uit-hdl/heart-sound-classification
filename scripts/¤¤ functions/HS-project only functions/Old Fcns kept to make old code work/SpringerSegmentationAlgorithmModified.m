function [deltaVec,qtVec,assigned_states] = SpringerSegmentationAlgorithmModified(id, IndAA, Fs, springerPar)

compDelta = false;
[~,~,~,heartRate,systolicTimeInterval] = findBestAuscArea(id,springerPar,Fs,compDelta);

deltaVec = zeros(1,4);
qtVec    = cell(1,4);
audioData = wav2TS(id,IndAA);
[PCG_Features, featuresFs] = getSpringerPCGFeatures(audioData, Fs);

for i=1:4
    [delta, ~, qt] = viterbiDecodePCG_Springer(PCG_Features, springerPar.piVector, ...
                    springerPar.Bmatrix,springerPar.totObsDist, heartRate(i), ...
                    systolicTimeInterval(i), featuresFs);
                
    deltaVec(i) = max(delta(end,:));
    qtVec{i}   = qt;
end

[~,bestAA] = max(deltaVec);
qt = qtVec{bestAA};
assigned_states = expand_qt(qt, featuresFs, Fs, length(audioData));

end