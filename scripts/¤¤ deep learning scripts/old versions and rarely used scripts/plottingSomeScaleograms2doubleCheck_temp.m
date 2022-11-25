data = HSdataTrainAndVal;
Nsamples = height(data);
Nds = 20;
% id = id_err(i_err);
k = randi(Nsamples);
id = HSdataTrainAndVal.UNIKT_LOPENR(k);
clf
for aa =1:4
    subplot(2,2,aa)
    x = wav2TS(id,aa);
    x = downsample(x,Nds);
    getScaleogram(x,1,true);
    hold on
    m = numel(x);
    assignedStates = getExpandedReprOfStates(nodes.loc{k,aa},nodes.state{k,aa},numel(x));
    plotAssignedStates(assignedStates,[2,4],['r','b'],1,40)
    title(sprintf('AR|MR|AS|MS = [%g,%g,%g,%g] murGrade=%g',...
    HSdata.ARGRADE_T72(HSdata.UNIKT_LOPENR==id),...
    HSdata.MRGRADE_T72(HSdata.UNIKT_LOPENR==id),...
    HSdata.ASGRADE_T72(HSdata.UNIKT_LOPENR==id),...
    HSdata.MSGRADE_T72(HSdata.UNIKT_LOPENR==id),...
    HSdata.(sprintf('Murmur_%g_grade_ref_ny_T72',aa))(HSdata.UNIKT_LOPENR==id)))
%     playHS(id,i)
%     pause(15)
end