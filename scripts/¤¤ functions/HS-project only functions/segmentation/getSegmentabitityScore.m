function segScore = getSegmentabitityScore(acf,acfInfo,NdsAcf)
% compute segmentation score based on the characteristics of the
% autocorellation function of the signal. acf is the autocorellation
% function, NdsAcf is how much to downsample the acf before estimating its
% timevarying mean, and acfInfo is a field that contains information about
% the peak used to identify the heart rate of the signal.

if nargin==2
    NdsAcf = 35;
end 

smoothFac = 3200/NdsAcf;
% fit only half of the acf:
x = downsample(acf(1:numel(acf)/2), NdsAcf);
% approximate the mean of the acf:
C = smoothdata(x,'gaussian',smoothFac);
% compute corrected for drifting mean:
z = x-C;
% get discrete cosine coefficients of the acf:
X = dct(z);
[~,ind] = sort(abs(X),'descend');
% approximate using the 4 most impactful coefficients:
NcosCoeff = 4;
X(ind(NcosCoeff+1:end)) = 0;
% take inverse transform to produce acf approximation:
xx = idct(X);
% find approximation error size:
RMSE = rms(xx-z);
RMSEbase = rms(z-C);
RMSEred = 100*(1-RMSE/RMSEbase);

% get the segmentation score
segScore = acfInfo.pp * RMSEred*(1.4*sigmoid(acfInfo.peakRatio - 1.4) + 1);

end