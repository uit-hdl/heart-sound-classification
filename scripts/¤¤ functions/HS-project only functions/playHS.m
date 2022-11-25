function player = playHS(id,aa)
x = audioread(sprintf('%.0f_hjertelyd_%g.wav',id,aa));
player = audioplayer(x,44100);
play(player);
end