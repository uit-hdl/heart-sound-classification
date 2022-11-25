function xCorrected = subtractTimeVarMean(x, smoothingFac)
% computes and subtracts slow changing time varying mean using gaussian
% smoothing.

if nargin==1
    smoothingFac = 3500;
end
% approximate slow time varying drift
C = smoothdata(x,'gaussian',smoothingFac);
xCorrected = x-C;
end
