function [x,fs] = wav2TS(id,auscLoc)
% reads the audio file and CONVERTS it to numeric vector. id is the
% LOPENUMMER.
[x,fs] = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,auscLoc));
end