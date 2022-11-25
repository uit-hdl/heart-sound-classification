function [heartRate,true_index,info] = findTrueHeartRate(acf, min_index, max_index, fs, thrHR, thrMult, plotIt)
% Estimates heart rate based on the auto correlation function acf, which is
% produced in the function "getHeartRateSchmidt". Fs is the sampling
% frequency (after downsampling, if it has been performed). min_index and
% max_index defines the ranges over which to search for the heart rate in
% the acf, and are obtained in "getHeartRateSchmidt". xc1 is the location
% of the largest peak, and xc2 is the location of the largest peak to the
% left of xc1, of which xc1 might be a multiple.

% In this function I offer an alternative way to estimate the heart rate.
% This method was created after I noticed that Schmidts method had a
% tendency to pick the wrong peak in the ACF to pick the wrong peak from
% which to calculate the heart rate, as the tallest peak is not always the
% one which corresponds to the heart rate. 
%% preliminary
if nargin==4
    thrHR   = .6;
    thrMult = .15;
    plotIt = false;
elseif nargin==5
    thrMult = .15;
    plotIt = false;
elseif nargin==6
    plotIt = false;
end

%%
% run below line to test code:
% acf = acFcn{1}; min_index = 0.5*fs; max_index = 2*fs; 

warning('off')
% C = smoothdata(acf,'gaussian',0.072563*Fs); % old code
xx = acf(min_index:max_index);
% xx = xx-C(min_index:max_index); % old code

% find three largest peaks over HR search range:
[pks,locs,~,p] = findpeaks(xx);
if isempty(pks)
    [yc1,xc1] = max(xx);
    pc1 = 0;
else
    % find largest peak:
    [yc1, i_max] = max(pks);
    xc1 = locs(i_max); % location of largest peak
    pc1 = p(i_max);    % peak prominence of largest peak
end

info.max_index = xc1 + min_index - 1;

%%% new way %%%
% set threshold above which a peak must exceed to qualify as candidate:
u_cand = thrHR*yc1;

% Find candidate peaks that are to the left of xc1 and that are above the
% peak consideration threshold:
indCand = and(locs < xc1, pks>u_cand);
cand_locs  = locs(indCand);
cand_pks   = pks(indCand);
cand_p     = p(indCand);

% get true (original coordinate system) locations of candidates and max
% peak:
cand_locsTrue = cand_locs + min_index - 1;
xc1           = xc1       + min_index - 1;

if isempty(cand_locs)
    % there are no peaks that satisfy conditions to be candidates; use xc1
    % to estimate heart rate:
    true_index = xc1;
    % ¤¤ gather analysis info ¤¤
    info.y_hr  = yc1;
    info.pp    = pc1; % peak prominence of HR-peak
    orderedPks = sort(pks);
    
    if numel(orderedPks)<=1
        info.peakRatio = nan;
    else
        info.peakRatio = yc1/orderedPks(end-1);
    end
    
else
    % there is atleast one candidate. find the candidate peak closest to
    % being a factor of xc1:
    x_thresholds = [1 - thrMult, 1 + thrMult]*(xc1/2);
    primeCandidates = isXinInterval(cand_locsTrue, x_thresholds);
    PrimeExists = sum(primeCandidates)~=0;
    
    
    if PrimeExists
        % in case more than one peak that satisfies prime criteria:
        [~,Ifactor] = max(cand_pks.*primeCandidates);
        xc2   = cand_locsTrue(Ifactor);
        yc2   = cand_pks(Ifactor);
        pc2   = cand_p(Ifactor);

        if plotIt
            figure
            plot(acf)
            hold on

            t = xc2;
            xx = xc1*.5;
            yy = yc1;
            line([xx - thrMult*t,xx + thrMult*t],[yy,yy]*0.6)
            line([xx - thrMult*t,xx - thrMult*t],[yy*0.6, 1])
            line([xx + thrMult*t,xx + thrMult*t],[yy*0.6, 1])
        end
        
        true_index = xc2;
        % ¤¤ gather analysis info ¤¤
        info.y_hr = yc2;
        info.pp    = pc2;
        info.peakRatio = yc2/yc1;
        
    else
        % factor test failed -- pick the largest peak in the fixed search
        % interval to estimate heart rate:
        true_index = xc1;
        % ¤¤ gather analysis info ¤¤
        info.y_hr = yc1;
        info.pp   = pc1;
        orderedPks = sort(pks);
        
        if numel(orderedPks)==1
                info.peakRatio = nan;
        else
            info.peakRatio = yc1/orderedPks(end-1);
        end
    end
end

% save some analytics
info.x_hr        = true_index/fs;
info.index_cand  = cand_locsTrue;
info.cand_pks    = cand_pks;
info.true_index  = true_index;
info.Pmax = [xc1/fs,yc1];

% compute heart rate
heartRate = 60/(true_index/fs);

warning('on')

end