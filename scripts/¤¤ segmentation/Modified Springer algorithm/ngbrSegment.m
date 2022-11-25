function [states,fitInfo] = ngbrSegment(id,Nds0,NdsAcf,HMMpar,segNgbr,segNoNgbr,plt,Fs)
% takes the 4 recordings corresponding to the identification number id and
% outputs segmentation of each using the confident neighbour method. Runs 
%%
if nargin==4
    segNgbr = true;
    segNoNgbr = false;
    plt.it = false;
    plt.compare = false;
    Fs = 44100;
elseif nargin==5
    segNoNgbr = false;
    plt.it = false;
    plt.compare = false;
    Fs = 44100;
elseif nargin==6
    plt.it = false;
    plt.compare = false;
    Fs = 44100;
elseif nargin==8
    Fs = 44100;
end

if isempty(Nds0)
    Nds0 = 20;
end
if isempty(NdsAcf)
    NdsAcf = 35;
end
if isempty(segNgbr)
    segNgbr = true;
end
if isempty(segNoNgbr)
    segNoNgbr = false;
end
%%
% get heart rate parameters estimates:
alfa = .22;
beta = .09;
[parEst,fitInfo] = getHRestUsingConfidentNgbr(id,alfa,beta,Nds0,NdsAcf);

if sum([segNgbr,segNoNgbr])==0
    % in case I accidently set both to false:
    segNgbr = true;
end
    
states = cell(1,4);
fs = floor(Fs/Nds0);
for aa=1:4
    x = downsample(wav2TS(id,aa),Nds0);
    
    if segNgbr
        % get segmentation using neighbourhood information:
        heartPar.rate        = parEst{1}(aa);
        heartPar.sysDuration = parEst{2}(aa);
        [assignedStatesNew] = runSpringerSegmentationAlgorithmMod(x, fs, HMMpar, heartPar);
        states{aa} = assignedStatesNew;
    end
    
    if segNoNgbr
        % get segmentation without using neighbourhood information:
        heartPar.rate        = fitInfo.h(aa);
        heartPar.sysDuration = fitInfo.s(aa);
    	[assignedStatesOld] = runSpringerSegmentationAlgorithmMod(x, fs, HMMpar, heartPar);
    end
    
    if plt.it 
        % plot the scaleogram with the states:
        subplot(2,2,aa)
        getScaleogram(x,1,true);
        hold on
        states2plot = [2,4];
        cols        = ['m','c'];
        plotAssignedStates(assignedStatesNew,states2plot,cols,1,48)
        title(sprintf('heartRate=%.3g, sys. dur. = %.2g',...
            heartPar.rate,heartPar.sysDuration))
        
        if plt.compare && segNoNgbr
            % plot the states computed by the old method:
            states2plot = [2,4];
            cols        = ['r','g'];
            plotAssignedStates(assignedStatesOld,states2plot,cols,1,50)
        end
    end
end

end