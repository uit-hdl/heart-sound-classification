% run this script first:
% plotAutoCorrelationFcn
%% FIND THE PEAK THAT GIVES HEART RATE
nDs = 35;
aa = 2;
x = downsample(acFcn{aa} - min(acFcn{aa}),nDs);
min_index = 0.5*fs/nDs;max_index = 2*fs/nDs;
C = smoothdata(x,'gaussian',3200/nDs);
clf
subplot(211)
    plot(x)
    hold on
    plot(C)
    xline([1/2,2]*fs/nDs)
    
subplot(212)
    x = x-C;
    [pks,locs,w,p] = findpeaks(x(min_index:max_index));
    [~,Ipeaks] = sort(p,'descend');
    Ilargest = Ipeaks(1:3);
    locs = locs(Ilargest);
    p = p(Ilargest);
    locs = locs + min_index - 1; % shift back to absolute (true) indeces
    pks  = pks(Ilargest);
    xc1 = locs(1);
    pc1 = p(1);
    plot(x)
    hold on
    xline([1/2,2]*fs/nDs)
    scatter(locs, pks)
    legend({'original signal','smoothed average'})
% hold on
% plot(C)
% subplot(212)
% plot(x-C)
% xline([1/2,2]*fs/nDs)
% legend({'signal minus timevarying mean'})

% we can use the peak prominence as a decision parameter for selecting
% which time series to use for segmentation. Greater peak prominence means
% better segementability.
% investigate largest peaks

% sort peaks and peak locations by locations (increasing order):
[~,ind] = sort(locs);
xxmax   = locs(ind);
pks     = pks(ind);
p       = p(ind);

f =@(t) abs(xxmax/t - round(xxmax/t));
% list of factors to try out
xvec = 20:.2:100; 
n = numel(xvec);
yvec = zeros(3,n);
% synergy contains values that measure the degree to which each factor
% could be a common factor.
synergy = zeros(2,n);
% CfScore shows how close to integer multiplicity the elements are
% for the given input
CfScore = zeros(1,n);
for i=1:n
    yvec(:,i) = f(xvec(i));
    % find synergy of the two best matches:
    [minVal,synergy(:,i)] = mink(yvec(:,i),2);
    CfScore(i) = max(minVal);
end

% clf
% plot(xvec,yvec)
% hold on 
% plot(xvec,CfScore,'k')

% VARIABLES
% Cf = common factor
% Rs = regularity score -- to what degree does it appear as a common factor


%%% the peaks that are closest to being integer multiples of common factor
% are: 

% regularity score: how close are the two elected peaks to being integer
% multiples of common factor.

% commonFactor: what is the common factor that represents the closest
% match.

% find out which Cf that gives best representation for the peak locations:
[Rs,I_Cf] = min(CfScore);
% get index
BestSynergy = synergy(:,I_Cf)
CfBest = xvec(I_Cf) %#ok<*NOPTS>
f(CfBest)


multiplicityScore = f(CfBest);
thrWeights = p.^2/mean(p.^2);
thrCf = 0.025;
testThr = thrCf * thrWeights
passTest = find(f(CfBest)<testThr)

IbestPeak = min(passTest);
bestPeak  = xxmax(IbestPeak);

if isempty(IbestPeak)
    overwiev.heartRate = 60/(xc1/fs)/nDs;
    overwiev.indBestPeak = ind(1);
    overwiev.bestPeak = xc1;
    overwiev.Cf       = CfBest;
    overwiev.thr      = testThr';
    overwiev.p        = p';
    overwiev.synScore = f(CfBest)'
else    
    overwiev.heartRate = 60/(bestPeak/fs)/nDs;
    overwiev.indBestPeak = IbestPeak;
    overwiev.bestPeak = bestPeak;
    overwiev.Cf       = CfBest;
    overwiev.p        = p';
    overwiev.thr      = testThr';
    overwiev.synScore = f(CfBest)'
end

