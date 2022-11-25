% In this script, I plot the scalograms and assigned states, as well as the
% autocorellation functions in 4 by 2 subplots. The acf's are plotted
% together with the discrete cosine functions, and the HR-peaks are marked
% with a dot whith a color which indicates the confidence score of the
% estimate.
%% preliminary
newk = 1;
k    = 0;
colMap = parula(255);
load HMMpar.mat
% select subset to investigate:
Jsubset = Jtrain0;

colormap parula
cmap = colormap;
%% Body
NdsAcf = 35;
Nds0 = 20;
fs0  = 44100;
fs = floor(fs0/Nds0);
redo  = 0;
kprev = k;

if newk==1
    k = k + 1 - redo; %#ok<*NASGU>
%     k = randi(height(data));
end
% good examples of cases where modified algorithm is improvement:
% 1746 1830 272 806 1623 1449 228 953
% examples of how things can go wrong:
% clear signal, but high irregularity: 89 265

id = HSdata.id(Jsubset(k));
Nds  = floor(30/Nds0);
% get true row index:
kk = find(findInd(HSdata.id,id));
% save heart parameters, and negative log likelihood
if redo==0
    saveHeart = zeros(4,2);
    NLL    = zeros(2,4);
    acFcn = cell(4,1);
end

% These are set to default:
% alfa = .22;
% beta = .11;
% Segment audio using the confident neighborhood method:
[states_ngbr,fitInfo] = ngbrSegment(id,[],[],HMMpar);
confidence = fitInfo.c;

close all
figure('units','normalized','outerposition',[0 0 1 1])
% *** plot scaleogram ***
for aa=1:4
    % read and downsample audio:
    [x,fs0] = wav2TS(id,aa);
    x = downsample(x,Nds0);
    fs = floor(fs0/Nds0);
    
    if redo==0
        heartPar.rate        = fitInfo.h(aa);
        heartPar.sysDuration = fitInfo.s(aa);
        [assignedStates,heartPar,acFcn{aa}] = runSpringerSegmentationAlgorithmMod(x, fs, HMMpar,heartPar);
    else
        [assignedStates,heartPar,acFcn{aa}] = runSpringerSegmentationAlgorithmMod(x, fs, HMMpar, best);
    end
    
    subplot(4,2,4+aa)
        scaleogramPlot(x,fs);
        hold on
        states2plot = [2,4];
        cols        = ['r','g'];
        murStr = sprintf(sprintf('Murmur_%g_grade_ref_ny_T72',aa));
        plotAssignedStates(assignedStates,states2plot,cols,fs,50)
        plotAssignedStates(states_ngbr{aa},states2plot,["magenta","cyan"],fs,100)

        if HSdata.(sprintf('noise%g',aa))(kk)==1
            murStr = '  noise';
        else
            murStr = sprintf('murGrade=%g', HSdata.(murStr)(kk) );
        end

        if aa==1
            s = sprintf('   (id=%s, k=%g)',string(id),kk);
        else
            s = '';
            ylabel ''
            xlabel ''
        end
        
        title(sprintf('HR=%.3g   systLen=%.2gs  %s%s',...
                    round(heartPar.rate),heartPar.sysDuration, murStr, s))

        saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
end

% *** plot ACF ***
plotSlowDrift = true;

for aa=1:4
    
    subplot(4,2,aa)
    % read and downsample audio:
    [x,fs0] = wav2TS(id,aa);
    x = downsample(x,Nds0);
    fs = floor(fs0/Nds0);
    
    [heartRate, systolicTimeInterval,acFcn{aa},acfInfo] = getHeartRateSpringerMod(x, fs);
    
    % choose new sampling rate for ACF before taking DCT approximation:
    NdsAcf = 35;
    fs_acf = round(fs/NdsAcf);
    acf = downsample(acFcn{aa}(1:floor(numel(acFcn{aa})/2)), NdsAcf);
    % get DCT approximation:
    n0 = 1
    [DCT, z, C] = DCTapproxOfACF(acf, fs, n0, 'smoothFac', 2*fs_acf);
    
    % find approximation error size:
    RMSE = rms(acf - (DCT.acf_reconst));
    RMSEbase = rms(acf-C);
    RMSEred = 100*(1-RMSE/RMSEbase);
    segScore = acfInfo.pp*RMSEred*(1.4*sigmoid(acfInfo.peakRatio - 1.4) + 1);
    

    % plot ACF:
    plot((1:numel(acf))/fs_acf, acf,'r','LineWidth',0.9)
    hold on
    % plot DCT reconstruction:
    plot((n0:numel(acf))/fs_acf, DCT.acf_reconst,'col',color2triplet("grey"))
    if plotSlowDrift
        plot((n0:numel(acf))/fs_acf, C,'blue')
    end
    
    % find color in colormap that corresonds to segmentation score:
    peakCol = cmap(floor(sigmoid((segScore-20)/30)*255),:);
    % plot max peak:
    scatter(acfInfo.Pmax(1), acfInfo.Pmax(2), [], peakCol,'filled')
    % plot identified peak:
    scatter(acfInfo.x_hr, acfInfo.y_hr, 100, peakCol,'Marker','p',...
                     'MarkerFaceColor',peakCol)
    % plot systole peak:
    scatter(acfInfo.sysLeng, acfInfo.y_syst, 30, peakCol,'Marker','o',...
                     'MarkerFaceColor','r')
                 
    % plot integer fraction search region:
    col = color2triplet("grey");
    xx = acfInfo.Pmax(1)/2;
    yy = acfInfo.Pmax(2);
    line([.85*xx,xx*1.15], [yy,yy]*0.6,'color',col)
    line([.85*xx,xx*.85],  [yy*0.6, 1],'color',col)
    line([1.15*xx,xx*1.15],[yy*0.6, 1],'color',col)
    
    % search interval lines
    xline(0.4,'col',color2triplet("grey"),'LineStyle','--')
    xline(2,'col',color2triplet("grey"),'LineStyle','--')
    
    % if the max peak gets replaced, plot following:
    if fitInfo.borrowStarEst(aa)==1
        a = scatter(acfInfo.Pmax(1), acfInfo.Pmax(2),100,'kx','LineWidth',1)
    end
    
    axis 'tight'
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
    coeffString = sprintf(' %g ',DCT.coeff);
    fitString   = sprintf(' %.2g ',fitInfo.z*100);
    xlim([0,4.0])
    
    title(sprintf('RMSEreduct=%.2g    peakProm=%.2g   peakRatio=%.2g  segScore=%.1f   \n borrowNgbrScore=%.2g   confidence=%g  ',...
        RMSEred, acfInfo.pp, acfInfo.peakRatio, segScore, fitInfo.z(aa),round(confidence(aa),1)))
end

fitInfo %#ok<*NOPTS>

