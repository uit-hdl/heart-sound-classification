load('Springer_B_matrix.mat');
load('Springer_pi_vector.mat');
load('Springer_total_obs_distribution.mat');
HMMpar.Bmatrix = Springer_B_matrix;
HMMpar.piVector = Springer_pi_vector;
HMMpar.totObsDist = Springer_total_obs_distribution;
% choose data frame:
% data = HSdataTrain;
% here you can select a smaller subset if desired:
Jsubset = 1:height(HSdata);
data = HSdata(Jsubset,:);
% data = HSdata(1:100,:);

%%
NdsAcf = 35; % downsampling of acf
Nds0 = 20;  % downsampling of signal
Fs = 44100; % sampling frequency
fs = floor(Fs/Nds0);
% define which set to explore;
set = data;
Nds = floor(30/Nds0);
n = height(data);

% ¤¤¤ CHOOSE SEGMENTATION ALGORITHM ¤¤¤
segAlg = "Springer";

stateInf.states = cell(n,1);
stateInf.id = zeros(n,1);
for i=1:n
    disp(i)
    id = set.id(i);
    
    if segAlg=="ngbrSegment"
        stateInf.states{i} = ngbrSegment(id,Nds0,NdsAcf,HMMpar); %#ok<*SAGROW>
        
    elseif segAlg=="Springer"
        for aa=1:4
            audio_data = wav2TS(id,aa);
            audio_data = downsample(audio_data,Nds0);
            stateInf.states{i}{aa} = runSpringerSegmentationAlgorithm(audio_data,...
                fs, HMMpar.Bmatrix, HMMpar.piVector, HMMpar.totObsDist); %#ok<*SAGROW>
        end
    end
    stateInf.id(i) = id;
end

%% Convert segmentation data into compact form
% clear nodes
cycles   = cell(n,4);
segLines = cell(n,4);
% nodesNew contains information about states in compact form
nodesNew.loc = cell(n,4);
nodesNew.state = cell(n,4);
nodesNew.seglines = cell(n,4);
% loc      = locations where the chain jumps to new state
% state    = which state chain jumped to
% seglines = positions where new cycle begins (end of diastole)
for i=1:n
    disp(i)
    for aa=1:4
        [cycles{i,aa},segLines{i,aa}] = ...
                            states2cycles(stateInf.states{i}{aa});
                        
        [nodesNew.loc{i,aa}, ...
         nodesNew.state{i,aa},...
         nodesNew.seglines{i,aa}] = ...
                     getCompactReprOfStates(stateInf.states{i}{aa});
    end
end
% save segmentation data:
% save('cyclesVal.mat','cyclesVal')
% save('segLinesVal.mat','segLinesVal')
% save('cyclesVal.mat','cyclesVal')

%% test
i = 4;
aa = 4;
Nds0 = 20;
id = stateInf.id(i);
x = downsample(wav2TS(id,aa),Nds0);
states = getExpandedReprOfStates(nodesNew.loc{i,aa},nodesNew.state{i,aa},numel(x));
close all
getScaleogram(x,[],true);
hold on
plotAssignedStates(states)
%% save as a field
% Run line to save (safety precaution to not overwrite existing file)
if 1==0
    save('nodes_Springer.mat','nodesNew')
end
% save('C:\Users\perwa\OneDrive - UiT Office 365\HSclassifProject\matlab files\saved variables\nodesAll.mat','nodesAll')
%% TEST
k = 3;
aa = 1;
id = set.id(k);
x = downsample(wav2TS(id,aa),Nds0);
S = ngbrSegment(id,Nds0,NdsAcf,HMMpar); %#ok<*SAGROW>
% [cyc,sL1] = states2cycles(states{k}{aa});
[cyc,sL2] = states2cycles(S{aa});

close all
M = getScaleogram(x,44100/Nds0,true);
hold on
plotAssignedStates(S{aa})
%%
load nodes

aa = 1;
close all
k = 1900;
figure
getScaleogramFromIndex(modelData,[],k,aa,true,nodes);

%%
load nodes
Jtraining = union(Jtrain0,Jval0);
Jtesting  = Jtest0;
nodes00.loc = cell(height(HSdata),4); 
nodes00.state = cell(height(HSdata),4); 
nodes00.seglines = cell(height(HSdata),4);

nodes00.loc(Jtraining,:) = nodes.loc;
nodes00.loc(Jtesting,:)  = nodesNew.loc;
nodes00.seglines(Jtraining,:) = nodes.seglines;
nodes00.seglines(Jtesting,:)  = nodesNew.seglines;
nodes00.state(Jtraining,:) = nodes.state;
nodes00.state(Jtesting,:)  = nodesNew.state;
%%
aa = 1;
close all
k = 1999;
figure
getScaleogramFromIndex(HSdata,[],k,aa,true,nodesAll);

