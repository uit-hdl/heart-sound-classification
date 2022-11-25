function [xx,F] = downSample(xx,F,Nds)
% function that takes cell arrays of vectors and downsamples them
xx{1} = downsample(xx{1},Nds);
xx{2} = downsample(xx{2},Nds);
if isvector(F{1})
    F{1}  = downsample(F{1} ,Nds);
    F{2}  = downsample(F{2} ,Nds);
else
    F{1}  = downsample(F{1}' ,Nds)';
    F{2}  = downsample(F{2}' ,Nds)';
end
end