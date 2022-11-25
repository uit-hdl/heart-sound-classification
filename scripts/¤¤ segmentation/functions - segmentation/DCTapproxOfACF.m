function [dctApprox,z,C] = DCTapproxOfACF(acf0,fs,n0,varargin)

%% set optional values
Ncoeff     = 4;
smoothFac  = 3200;
plotApprox = false;
%% optional name-value-pair arguments:
p = inputParser;
addOptional(p, 'Ncoeff'    , Ncoeff    , @(x) isnumeric(x));
addOptional(p, 'smoothFac' , smoothFac , @(x) isnumeric(x));
addOptional(p, 'plotApprox', plotApprox, @(x) islogical(x));
parse(p,varargin{:});
% detach:
Ncoeff     = p.Results.Ncoeff;
smoothFac  = p.Results.smoothFac;
plotApprox = p.Results.plotApprox;

%% *** Discrete Cosine Transform ***

% *** preprocessing the ACF ***
% Remove first n0 elements of acf:
acf = acf0(n0:end);
% approximate the slow drift in mean:
C = smoothdata(acf,'gaussian',smoothFac);
% compute ACF corrected for drifting mean which is to be approximated by
% DCT:
z = acf-C;

% *** approximate the ACF using the DCT ***
% get the DCT coefficients of the acf:
dct_coeff = dct(z);
N = numel(z);
% sort the coefficients by the magnitudes:
[dct_coeff_sorted,J_sorted] = sort(abs(dct_coeff),'descend');
% set the other coefficients to zero:
dct_coeff(J_sorted(Ncoeff+1:end)) = 0;
% get index of most impactful coefficients:
J_approx = J_sorted(1:Ncoeff);
% slowestFreq = sum(J_approx .* XX(1:NcosCoeff)/sum(XX(1:NcosCoeff)));
% take inverse	transform to produce acf approximation:
z_reconst = idct(dct_coeff);

% compute the frequency corresponding to each coefficient:
fn = pi/(2*N)*(2*J_approx - 1);
% Number of time increments to finish one cycle whe frequency is fn:
Tn = 2*pi./fn;

% *** convert period to seconds and frequency to 1/s:
T_dct = Tn'/fs;
f_dct = 1./T_dct;


% ### collect output in a structure ###
dctApprox.f     = f_dct;
dctApprox.T     = T_dct;
dctApprox.coeff = dct_coeff_sorted(1:Ncoeff)';
dctApprox.acf_reconst = z_reconst + C;

% *** plot:
if plotApprox
    plot((1:numel(z))/fs, acf)
    hold on
    plot((1:numel(z_reconst))/fs, z_reconst)
end


end

% find approximation error size:
% RMSE = rms(xx-z);
% RMSEbase = rms(z-C);
% RMSEred = 100*(1-RMSE/RMSEbase);
% segScore = acfInfo.pp*RMSEred*(1.4*sigmoid(acfInfo.peakRatio - 1.4) + 1);
% plot((1:numel(z))*NdsAcf,idct(X)+C,'r')
% hold on
% % find color in colormap that corresonds to segmentation score:
% peakCol = colMap(floor(sigmoid((segScore-15)/30)*255),:);
% scatter(acfInfo.true_index, acf(round(acfInfo.true_index/NdsAcf)),[],...
%                                         peakCol,'filled')
% if fitInfo.borrowStarEst(aa)==1
% scatter(acfInfo.true_index, acf(round(acfInfo.true_index/NdsAcf)),70,'rx')
% end
% axis 'tight'
% saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
% coeffString = sprintf(' %g ',sort(approxCoeff));
% %     fitString = sprintf(' %.2g ',fitInfo.z*100);
% title(sprintf('RMSEreduction=%.2g   peakProminence=%.2g   peakRatio=%.2g   segScore=%.2g   \n borrowNgbrScore=%.2g  ',...
% RMSEred, acfInfo.pp, acfInfo.peakRatio, segScore, fitInfo.z(aa)))