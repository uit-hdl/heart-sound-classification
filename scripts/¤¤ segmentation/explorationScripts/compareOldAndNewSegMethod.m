% In this script I compare the springer segmentation algorithm to my
% modified version (modification 1, individual segmentation).

% load indeces for which disagreement occurs:
load('kComp.mat');

data = HSdataTrain;
newk = 1;
k=1;
%% PRELIMINARY
Nds0 = 20;
fs = floor(44100/Nds0);
redo = 0;
if newk==1
    k = k + 1 - redo;
end
id = data.UNIKT_LOPENR(kComp(k)); %#ok<*NASGU>
Nds  = floor(30/Nds0);

% save heart parameters, and negative log likelihood
if redo==0
    saveHeart = zeros(4,2);
    NLL    = zeros(2,4);
    acFcn = cell(4,1);
end
%% COMPARE SCALEOGRAMS
if newk==1
    k = k + 1 - redo;
end
id = data.UNIKT_LOPENR(kComp(k));
Nds  = floor(30/Nds0);
clf
%%% New Method %%%
for aa=1:4
    x = downsample(wav2TS(id,aa),Nds0);
    if redo==0
        [assignedStates,heartPar] = runSpringerSegmentationAlgorithmMod(x, fs, HMMpar);
    else
        [assignedStates,heartPar] = runSpringerSegmentationAlgorithmMod(x, fs, HMMpar,best);
    end
    x = downsample(x,Nds);
    
    subplot(4,2,1 + 2*(aa-1))
    getScaleogram(x,1,true);
    hold on
    states2plot = [2,4];
    cols        = ['r','g'];
    plotAssignedStates(assignedStates,states2plot,cols,Nds,45)
    if aa==1
        title('modified method')
    end
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
end

%%% old method %%%
for aa=1:4
    x = downsample(wav2TS(id,aa),Nds0);
    if redo==0
        [assignedStates,heartPar,NLL(1,aa)] = runSpringerSegmentationAlgorithmOrgl(x, fs, HMMpar);
    else
        [assignedStates,heartPar,NLL(2,aa)] = runSpringerSegmentationAlgorithmOrgl(x, fs, HMMpar,best);
    end
    x = downsample(x,Nds);
    
    subplot(4,2,2*aa)
    getScaleogram(x,1,true);
    hold on
    states2plot = [2,4];
    cols        = ['r','g'];
    plotAssignedStates(assignedStates,states2plot,cols,Nds,45)
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
    if aa==1
        title('old method')
    end
end

%% COMPARE 1 SEGMENTATION ATTEMPT AT A TIME
if newk==0
    k = k + 1;
%     kk = kComp(k);
    kk=k;
end
% kk = 146 is a clear example of the bug that the modified version fixes
kk = 146;
id = data.UNIKT_LOPENR(kk);
Nds  = floor(30/Nds0);
clf
%%% New Method %%%
for i=1:2
    x = downsample(wav2TS(id,1),Nds0);
    if i==1
        [assignedStates,heartPar] = runSpringerSegmentationAlgorithmMod(x, fs, HMMpar);
    else
        [assignedStates,heartPar] = runSpringerSegmentationAlgorithmOrgl(x, fs, HMMpar);
    end
    x = downsample(x,Nds);
    
    subplot(2,2,i)
    getScaleogram(x,1,true);
    hold on
    states2plot = [2,4];
    cols        = ['r','g'];
    plotAssignedStates(assignedStates,states2plot,cols,Nds,45)
    title(sprintf('method %g',3-i))
end

for i=1:2
    subplot(2,2,2+i)
    x = downsample(wav2TS(id,1),Nds0);
    x = downsample(x,Nds);
    if i==1
        [heartRate, systolicTimeInterval,acFcn{i},p] = getHeartRateSpringerMod(x, fs, true);
    else
        [heartRate, systolicTimeInterval,acFcn{i}] = getHeartRateSchmidtOrgl(x, fs, true);
    end
    hold on
    
    NdsAcf = 35;
    smoothFac = 3200/NdsAcf;    
    x = downsample(acFcn{i}(1:numel(acFcn{i})/2),NdsAcf);
    C = smoothdata(x,'gaussian',smoothFac);
    % compute corrected for drifting mean:
    x = x-C;
    % get discrete cosine coefficients of the acf:
    X = dct(x);
    [XX,ind] = sort(abs(X),'descend');
    % approximate using the 4 most impactful coefficients:
    NcosCoeff = 4;
    X(ind(NcosCoeff+1:end)) = 0;
    approxCoeff = ind(1:NcosCoeff);
    slowestFreq = sum(approxCoeff .* XX(1:NcosCoeff)/sum(XX(1:NcosCoeff)));
    % take inverse	ransform to produce acf approximation:
    xx = idct(X);
    % find approximation error size:
    RMSE = rms(xx-x);
    RMSEbase = rms(x-C);
    RMSEred = 100*(1-RMSE/RMSEbase);
%     plot((1:numel(x))*35,idct(X)+C,'r')
    axis 'tight'
    saveHeart(i,:) = [heartPar.rate,heartPar.sysDuration];
    coeffString = sprintf(' %g ',sort(approxCoeff));
end

%% COMPARE ALL

clf
%%% New Method %%%
for aa=1:4
    x = downsample(wav2TS(id,aa),Nds0);
    if redo==0
        [assignedStates,heartPar,NLL(1,aa)] = runSpringerSegmentationAlgorithm(x, fs, HMMpar);
    else
        [assignedStates,heartPar,NLL(2,aa)] = runSpringerSegmentationAlgorithm(x, fs, HMMpar,best);
    end
    x = downsample(x,Nds);
    
    subplot(4,4,4+aa)
    getScaleogram(x,1,true);
    hold on
    states2plot = [2,4];
    cols        = ['r','g'];
    plotAssignedStates(assignedStates,states2plot,cols,Nds,50)
    
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
end

for aa=1:4
    subplot(4,4,aa)
    x = downsample(wav2TS(id,aa),Nds0);
    x = downsample(x,Nds);
    [heartRate, systolicTimeInterval,acFcn{aa},p] = getHeartRateSchmidtV2(x, fs, true);
    hold on
    
    NdsAcf = 35;
    smoothFac = 3200/NdsAcf;    
    x = downsample(acFcn{aa}(1:numel(acFcn{aa})/2),NdsAcf);
    C = smoothdata(x,'gaussian',smoothFac);
    % compute corrected for drifting mean:
    x = x-C;
    % get discrete cosine coefficients of the acf:
    X = dct(x);
    [XX,ind] = sort(abs(X),'descend');
    % approximate using the 4 most impactful coefficients:
    NcosCoeff = 4;
    X(ind(NcosCoeff+1:end)) = 0;
    approxCoeff = ind(1:NcosCoeff);
    slowestFreq = sum(approxCoeff .* XX(1:NcosCoeff)/sum(XX(1:NcosCoeff)));
    % take inverse	ransform to produce acf approximation:
    xx = idct(X);
    % find approximation error size:
    RMSE = rms(xx-x);
    RMSEbase = rms(x-C);
    RMSEred = 100*(1-RMSE/RMSEbase);
    plot((1:numel(x))*35,idct(X)+C,'r')
    axis 'tight'
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
    coeffString = sprintf(' %g ',sort(approxCoeff));
end

%%% old method %%%

for aa=1:4
    x = downsample(wav2TS(id,aa),Nds0);
    if redo==0
        [assignedStates,heartPar,NLL(1,aa)] = runSpringerSegmentationAlgorithmOrgl(x, fs, HMMpar);
    else
        [assignedStates,heartPar,NLL(2,aa)] = runSpringerSegmentationAlgorithmOrgl(x, fs, HMMpar,best);
    end
    x = downsample(x,Nds);
    
    subplot(4,4,8+4+aa)
    getScaleogram(x,1,true);
    hold on
    states2plot = [2,4];
    cols        = ['r','g'];
    plotAssignedStates(assignedStates,states2plot,cols,Nds,50)
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
end

for aa=1:4
    subplot(4,4,8+aa)
    x = downsample(wav2TS(id,aa),Nds0);
    x = downsample(x,Nds);
    [heartRate, systolicTimeInterval,acFcn{aa}] = getHeartRateSchmidtOrgl(x, fs, true);
    hold on
    
    NdsAcf = 35;
    smoothFac = 3200/NdsAcf;    
    x = downsample(acFcn{aa}(1:numel(acFcn{aa})/2),NdsAcf);
    C = smoothdata(x,'gaussian',smoothFac);
    % compute corrected for drifting mean:
    x = x-C;
    % get discrete cosine coefficients of the acf:
    X = dct(x);
    [XX,ind] = sort(abs(X),'descend');
    % approximate using the 4 most impactful coefficients:
    NcosCoeff = 4;
    X(ind(NcosCoeff+1:end)) = 0;
    approxCoeff = ind(1:NcosCoeff);
    slowestFreq = sum(approxCoeff .* XX(1:NcosCoeff)/sum(XX(1:NcosCoeff)));
    % take inverse	ransform to produce acf approximation:
    xx = idct(X);
    % find approximation error size:
    RMSE = rms(xx-x);
    RMSEbase = rms(x-C);
    RMSEred = 100*(1-RMSE/RMSEbase);
    plot((1:numel(x))*35,idct(X)+C,'r')
    axis 'tight'
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
    coeffString = sprintf(' %g ',sort(approxCoeff));
end
