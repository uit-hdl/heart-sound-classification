
subplot(211)
[cD cA] = getDWT(x,200,'morl');
imagesc(abs(cD))
subplot(212)
getScaleogram(x,1,true);
