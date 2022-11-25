clf
xx = 1:1:255; 
for i=1:numel(xx)
plot(xx(i),1,'*','color',[xx(i),0,0]/255);
hold on
end

%%
clf
colormap parula 
colMap = parula(255)
for i=1:250
    scatter(1,2,[],colMap(i,:),'filled')
    pause(.01)
end

clf
N = 10;
x = linspace(0,3*pi,N);
y = cos(x) + rand(1,N);
c = linspace(1,10,length(x));
scatter(x,y,[],c)

%%
clf
x = downsample(wav2TS(id,aa),Nds0);
    
    x = downsample(x,Nds);
    [heartRate, systolicTimeInterval,acFcn{aa},info] = getHeartRateSchmidt(x, fs, true);
    p = info.pp;
    hold on
    
    NdsAcf = 1;
    smoothFac = 3200/NdsAcf;    
    x = downsample(acFcn{aa}(1:floor(numel(acFcn{aa})/2)),NdsAcf);
    x = acFcn{aa}(1:floor(numel(acFcn{aa})/2));
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
    segScore = p*RMSEred;
    plot((1:numel(x))*NdsAcf,idct(X)+C,'r')
    plot(x)
    hold on
     % find color in colormap that corresonds to segmentation score:
    peakCol = colMap(floor(sigmoid((segScore-15)/30)*255),:);
    scatter(info.true_index, x(floor(info.true_index))+C,[],...
                                                peakCol,'filled')
    if fitInfo.borrowStarEst(aa)==1
        scatter(info.true_index, x(info.true_index)+C,70,'rx')
    end
    axis 'tight'
    saveHeart(aa,:) = [heartPar.rate,heartPar.sysDuration];
    coeffString = sprintf(' %g ',sort(approxCoeff));
%     fitString = sprintf(' %.2g ',fitInfo.z*100);
    title(sprintf('RMSEreduction=%.2g   peakProminence=%.2g   peakRatio=%.2g   segScore=%.2g   \n borrowNgbrScore=%.2g  ',...
        RMSEred, p, info.peakRatio, segScore, fitInfo.z(aa)))

