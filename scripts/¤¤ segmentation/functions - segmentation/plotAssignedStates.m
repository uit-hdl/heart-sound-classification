function plotAssignedStates(assignedStates,states2plot,colors,fs,level)
% function that plots the assigned states (expanded form) specified by
% states2plot using the colors specified in the string vector cols. Nds is
% the downsampling factor of the data. Assumes that you have already
% plotted the scaleogram of the data. If Nds==20 (the default), set Nds to
% 1. This does not make sense, but it works, and I am too lazy to fix this
% bug every single script where it is used...

%% *** example input ***
% assignedStates = [1,1,1,1,2,2,2,2,3,3,3,3,4,4,4,1,1,1,1...]
% states2plot    = [2,4];
%%
if nargin==1
    states2plot = [2,4];
    colors  = ['r','b'];
    fs    = 1;
    level = 40;
elseif nargin==2
    colors  = ['r','b'];
    fs   = 1;
    level = 40;
elseif nargin==3
    fs   = 1;
    level = 40;
elseif nargin==4
    level = 40;
end

if isempty(states2plot)
    % by default plots systole and diastole:
    states2plot = [2,4];
end
if isempty(colors)
    colors = ['r','b'];
end

%% function body

for i=1:numel(states2plot)
    Jstate = find(assignedStates==states2plot(i));
    scatter(Jstate/fs, ones(1,length(Jstate))*level,colors(i),'.');
end

end