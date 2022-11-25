%%% script where I import data, and introduce convenience variables to make
%%% the data easier to work with.
load('seta')
load('setb')
load('seta_timing')
             
%%% S1 and S2 timing data
loc   = seta_timing.location;
sound = seta_timing.sound;

% names of sound files that have labelled S1 and S2
namesLabeledSignals = unique(seta_timing.fname);
% number of files with labeled peaks
Na.lab = length(namesLabeledSignals);
Ia.lab = zeros(Na.lab,length(seta.fname));
for k=1:Na.lab
    Ia.lab(k,:) = namesLabeledSignals(k)==seta.fname;
end

% Indeces of signals for which the peaks have been labelled
Ia.lab = sum(Ia.lab);     Ja.lab = find(Ia.lab);

%%% store the locations of the labelled peaks in cells
Ja.S  = cell(Na.lab,1);
Ja.S1 = cell(Na.lab,1);
Ja.S2 = cell(Na.lab,1);

for k=1:Na.lab
    I = seta.fname(Ja.lab(k))==seta_timing.fname;
    I_S1 = and(I,sound=="S1"); % indeces for S1 sounds
    I_S2 = and(I,sound=="S2");
    
    Ja.S{k}  = loc(I);
    Ja.S1{k} = loc(I_S1);
    Ja.S2{k} = loc(I_S2);
end
clear('loc','sound','I','I_S1','I_S2','k')

%%% INDEX VECTORS  
%     ** LOGICAL **                    ** REGULAR **
Ia.nor = seta.label=="normal";        Ja.nor      = find(Ia.nor);
Ia.mur = seta.label=="murmur";        Ja.mur      = find(Ia.mur);
Ia.art = seta.label=="artifact";      Ja.art      = find(Ia.art);
Ia.exthls = seta.label=="extrahls";   Ja.exthls   = find(Ia.exthls);
Ia.train = ~(seta.label == "");       Ja.train    = find(Ia.train);

Ib.nor = setb.label=="normal";           Jb.nor        = find(Ib.nor);
Ib.mur = setb.label=="murmur";           Jb.mur        = find(Ib.mur);
Ib.extrast = setb.label=="extrastole";   Jb.extrast    = find(Ib.extrast);
Ib.train = ~(setb.label == "");          Jb.train      = find(Ib.train);
Ib.noisy = or(setb.sublabel=="noisynormal"...
             ,setb.sublabel=="noisymurmur"); Jb.noisy  = find(Ib.noisy);
Ib.norOrMur = or(Ib.nor,Ib.mur);         Jb.norOrMur   = find(Ib.norOrMur);

%%% AMOUNTS OF EACH DATATYPE
Na.full = length(seta.fname); % total amount of data type (A)
Nb.full = length(setb.fname); % total amount of data type (A)
Na.train = sum(Ia.train); 
Nb.train = sum(Ib.train); 
Na.nor = sum(Ia.nor);
Nb.nor = sum(Ib.nor);
Na.mur = sum(Ia.mur);
Nb.mur = sum(Ib.mur);
Nb.norOrMur = sum(Ib.norOrMur);
Na.art = sum(Ia.art);
Na.exthls = sum(Ia.exthls);
Nb.extrast = sum(Ib.extrast);
Nb.noisy = sum(Ib.noisy);

%%% CREATE CELLS WHICH CONTAINS ALL THE SIGNALS

Xa = cell(Na.full,1);
Xb = cell(Nb.full,1);

% CONVERT INTPUT DATA FROM WAV. FORMAT TO TIMESERIES:
for k=1:Na.full
    Xa{k,1} = audioread(seta.fname(k));
end
for k=1:Nb.full
    Xb{k,1} = audioread(setb.fname(k));
end



