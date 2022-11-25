%% Guess case (training)
toggle=0;
toggle = toggle + 1;
toggle = mod(toggle, 2);
if toggle==1
    r = rand(1);
    if r>.5
        Ind = find(stats.I.clinSignVHD.any);
    else
        Ind = find(~stats.I.clinSignVHD.any);
    end
    k = randsample(height(Ind),1);
    k = Ind(k);
end
k=0;
k=k+1;
% choose subset to investigate
% Ind = stats.I.MS.GEQ.grade1.any;
Ind = and(~HSdata.murGradeSum>0,HSdata.ARgrade==4);
% Ind = stats.I.innocentMurm.any;
set = HSdata(Ind,:);
id = set.UNIKT_LOPENR(k);
% murmur grade
I = plotValCVTrain.UNIKT_LOPENR==id;

g(1) = plotValCVTrain.ARGRADE_T72(I);
g(2) = plotValCVTrain.MRGRADE_T72(I);
g(3) = plotValCVTrain.ASGRADE_T72(I);
g(4) = plotValCVTrain.MSGRADE_T72(I);

% figure()
% for i=1:4
% subplot(2,2,i)
%     plot(downsample(x,1+(i-1)*20))
% end


% clf
% for i=1:4
% x = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,1));
% % x = downsample(x,20);
% subplot(2,2,i)
%     plot(lowpass(downsample(x,1+(i-1)*10),0.25))
% end
% clf
% for i=1:4
% x = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,1));
% x = downsample(x,20);
% subplot(2,2,i)
%     plot(lowpass(x,1/(2*i+1)))
% end

clf
M = cell(1,4);
for i=1:4
str = sprintf('%g-', g);
x = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,i));
% player = audioplayer(x,44100);
% play(player)
x = downsample(x,20);
x = lowpass(x,0.1);
Fs = 44100*20;
subplot(2,2,i)
    [cfs,~] = cwt(x,'amor',Fs);
    tms = (0:length(x)-1)/Fs; % sampling times
    M{i} = abs(cfs);
    M{i} = M{i}(1:60,:);
%     M{i} = pseudoLog(M{i},1,inf);
    image(M{i}/0.00015)
    
    title(sprintf('age%g, PO2=%.2g, dyspneaFastUphill=%.2g, BMI=%.2g, \n std=%.2g, avgSize=%.2g',set.AGE_T7(k),...
        set.PO2_T72(k)/mean(set.PO2_T72,'omitnan'),...
        set.DYSPNEA_FAST_UPHILL_T7(k),...
        set.BMI_T7(k)/mean(set.BMI_T7,'omitnan'),...
        std(x), mean(abs(x)) ))
    
    if toggle==0
        title(sprintf('AR|MR|AS|MS: %s \n murgrade%g',str,mg(i)))
    
    if i==1
        title(sprintf('AR|MR|AS|MS: %s \n murgrade%g',str,mg(i)))
    elseif i==2
        title(sprintf('BMI: %.2g \n murgrade%g',plotValCVTrain.BMI_T7(I),mg(i)))
    end
    end
end
x = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,1));
player = audioplayer(x,44100);
play(player)


%%
[m,ind] = max(M{1});
clf
plot(m)

%% image analysis; a first attempt
%%% smoothen and downsample image %%%
t = 2;
Mds  = downsample(M{t}' ,20)'  ;
Mds  = imgaussfilt(Mds,3);
sz = size(Mds);
clf
subplot(311)
    imagesc(M{t})
    axis 'tight'
subplot(312)
    imagesc(Mds)
    axis 'tight'
subplot(313)
    [m,ind] = max(Mds);
    plot(m)
    peaks = findpeaks(m);
    axis 'tight'
%%
xx = reshape(Mds, [sz(1)*sz(2),1]); % reshape image columnwise into a vector

K = 3; % number of classes
[theta,prior,p,samples] = normmix_gibbs(xx,K,[],0);
figure
for i=1:K
    subplot(3,2,i)
    p_im = reshape(p(:,i),[sz(1),sz(2),1]);
    imagesc(rgbimage(p_im))
    title(sprintf('certainty class %d', i))
end
[cl,cl_ind,p] = normmix_classify(xx,theta,prior);

class_im = reshape(cl_ind,[sz(1),sz(2),K]);
subplot(324)
image(rgbimage(class_im))
title(sprintf( "classification with %g classes",K))
subplot(326)
imagesc(Mds)
title(sprintf( "classification with %g classes",K))
%%
% subplot(337)
% imagesc(rgbimage(x(:,:,[2 3 4])))
% title "original image"
