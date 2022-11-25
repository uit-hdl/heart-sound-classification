% In this script I plot scaleograms along with assigned states in a 2 by 2
% plot. Included in the subtitles are noise indicators and murmur-grade.
data = HSdata;

noiseLvl{1} = data.NOISE_NUMBER_1234_REF_T72;
noiseLvl{2} = data.MURMUR_2NOISE_REF_T72;
noiseLvl{3} = data.MURMUR_3NOISE_REF_T72;
noiseLvl{4} = data.MURMUR_4NOISE_REF_T72;

clearvars I
% agreed noise in atleast one location
I.noiseAny = (noiseLvl{1}+noiseLvl{2}+noiseLvl{3}+noiseLvl{4})>0;

%% Get PCG Features:
if ~exist('k','var')
    k = 0;
end

Nds0 = 20;
fs = floor(44100/Nds0);
% define which set to explore;
set = data;
redo = 0;
k = k + 1 - redo;
id = set.UNIKT_LOPENR(k);
Nds  = floor(30/Nds0);
clf
saveHeart = zeros(4,2);

for aa=1:4
    x = downsample(wav2TS(id,aa),Nds0);
    if redo==0
        [assignedStates,heartPar] = runSpringerSegmentationAlgorithmOrgl(x, fs, HMMpar);
    else
        [assignedStates,heartPar] = runSpringerSegmentationAlgorithm(x, fs, HMMpar,best);
    end
    x = downsample(x,Nds);
    subplot(2,2,aa)
    % index AuscultationArea

    getScaleogram(x,1,true);
    hold on
    states2plot = [2,4];
    cols        = ['r','g'];
    plotAssignedStates(assignedStates,states2plot,cols,Nds,40)
    title(sprintf('heartRate=%.2g, sys. dur. = %.2g, noiseLvl=%.1g, murGrade=%g',...
        heartPar.rate,heartPar.sysDuration,noiseLvl{aa}(k),data.maxMeanMurGrade(k)))
    
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
end
aa   = 1;
best.rate        = saveHeart(aa,1);
best.sysDuration = saveHeart(aa,2);

% It appears that it is not a good idea to assume that heart rate is
% constant between recordings. For k = 15, we get a heart rate of 120 in the
% mitral position, yet in the other 3 positions we get a heart rate of 80.
