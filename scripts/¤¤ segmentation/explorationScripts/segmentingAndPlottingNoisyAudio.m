loadHMMpar

aa = 4;
k = 1;
Inoise = and(HSdata.(sprintf('noiseOnly',aa)), HSdata.ASgrade>0); %#ok<*CTPCT>
Jnoise = find(Inoise);
k0 = Jnoise(k);
id = HSdata.id(k0);

[x,Fs0] = wav2TS(id,aa);
Nds = 20;
x = downsample(x,20);
Fs = floor(Fs0/Nds);

id2ind(id,HSdata)
[assignedStates, heartPar, acf, info] = runSpringerSegmentationAlgorithmMod...
                                           (x, Fs, HMMpar, [], true)

close all
subplot(211)
    getScaleogram(x,[],true);
    hold on
    plotAssignedStates(assignedStates)
    title(sprintf('AR=%g,MR=%g,AS=%g,MS=%g',HSdata.ARgrade(k0),HSdata.MRgrade(k0),...
                                            HSdata.ASgrade(k0),HSdata.MSgrade(k0)) )
subplot(212)
    plot(acf)


