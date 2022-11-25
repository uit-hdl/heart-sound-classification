% k=0
% k=k+1;

% choose subset to investigate
% Ind = stats.I.MR.GEQ.grade3.any;
% Ind = and(~stats.I.audMurWeak.any, stats.I.AR.GEQ.grade4.any);
Ind = HSdata.MRGRADE_T72>=4;

% Ind = stats.I.innocentMurm.any;
data = HSdata(Ind,:);
id = data.UNIKT_LOPENR(k);
% id = set.UNIKT_LOPENR(k);
% murmur grade
mg = HSdataTrain{Ind,258:261};
I  = HSdataTrain.UNIKT_LOPENR==id;

g(1) = HSdataTrain.ARGRADE_T72(I);
g(2) = HSdataTrain.MRGRADE_T72(I);
g(3) = HSdataTrain.ASGRADE_T72(I);
g(4) = HSdataTrain.MSGRADE_T72(I);

% figure()
% for i=1:4
% subplot(2,2,i)
%     plot(downsample(x,1+(i-1)*20))
% end


clf
for i=1:4
x = wav2PCG(id,i);
x = downsample(x,20);
subplot(2,2,i)
    plot(downsample(x,1+(i-1)*20))
end

clf
M = cell(1,4);
for i=1:4
str = sprintf('%g-',g);
x = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,i));
% player = audioplayer(x,44100);
% play(player)
x = downsample(x,20);
Fs = 44100*40;
subplot(2,2,i)
    [cfs,~] = cwt(x,'amor',Fs);
    tms = (0:length(x)-1)/Fs; % sampling times
    M{i} = abs(cfs);
%     M{i} = pseudoLog(M{i},1,inf);
    imagesc(M{i}/0.00015)
    title(sprintf('age%g, PO2=%.2g, dyspneaFastUphill=%g, BMI=%.2g',set.AGE_T7(k),...
        set.PO2_T72(k)/mean(set.PO2_T72,'omitnan'),...
        set.DYSPNEA_FAST_UPHILL_T7(k),...
        set.BMI_T7(k)/mean(set.BMI_T7,'omitnan')))

        title(sprintf('AR|MR|AS|MS: %s \n murgrade%g',str,mg(i)))
    
    if i==1
        title(sprintf('AR|MR|AS|MS: %s \n murgrade%g',str,mg(i)))
    elseif i==2
        title(sprintf('BMI: %.2g \n murgrade%g',HSdataTrain.BMI_T7(I),mg(i)))
    end

end
x = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,2));
player = audioplayer(x,44100);
play(player)