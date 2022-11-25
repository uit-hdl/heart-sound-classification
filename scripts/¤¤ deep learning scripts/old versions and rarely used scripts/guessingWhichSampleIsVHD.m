%% description
% in this script I try to train my ability to guess whether a recording
% comes from a sick or healthy person.
%% initialize
ticker = 0;
%%
clearvars id N k x I J
%%
load nodesAll
data = HSdata;
I.sick    = data.MRGRADE_T72==4;
J.sick    = find(I.sick);
I.control = data.MRGRADE_T72==0;
J.control = find(I.control);
N.sick = sum(I.sick);
N.control = sum(I.control);
k.sick = 1;
k.control = randi(N.control);
ticker = 0;
t = randperm(2);
%%
aa  = 4;
Nds = 20;
ticker = mod(ticker+1,2);

k.sick = k.sick*(ticker==0) + (k.sick + 1)*(ticker==1);
kk.sick = J.sick(k.sick);
k.control = k.control*(ticker==0) + randi(N.control)*(ticker==1);
kk.control = J.control(k.control);

id.sick = data.UNIKT_LOPENR(kk.sick);
id.control = data.UNIKT_LOPENR(kk.control);
x.sick = wav2TS(id.sick,aa);
x.sick = downsample(x.sick,Nds);
x.control = wav2TS(id.control,aa);
x.control = downsample(x.control,Nds);

gradeAS = data.ASGRADE_T72(kk.sick);
gradeMS = data.MSGRADE_T72(kk.sick);

states.sick = getExpandedReprOfStates(...
    nodes.loc{kk.sick,aa}, nodes.state{kk.sick,aa}, numel(x.sick));
states.control = getExpandedReprOfStates(...
    nodes.loc{kk.control,aa}, nodes.state{kk.control,aa}, numel(x.control));

labels = {'sick','not sick'};
t = (ticker==0)*t + (ticker==1)*randperm(2);
if ticker==1
    clf
end
subplot(2,1,t(1))
    getScaleogram(x.sick,1,true);
    hold on
    plotAssignedStates(states.sick,[2,4],['r','b'],1)
    if ticker==0
        title(sprintf('status: %s,k=%g,ASgrade=%g,MSgrade=%g','sick',kk.sick,gradeAS,gradeMS));
    end
subplot(2,1,t(2))
    getScaleogram(x.control,1,true);
    hold on
    plotAssignedStates(states.control,[2,4],['r','b'],1)
if ticker==0
    title(sprintf('status: %s','healthy'));
end

%%
k=1
id = data.UNIKT_LOPENR(k);
x = wav2TS(id,aa);





