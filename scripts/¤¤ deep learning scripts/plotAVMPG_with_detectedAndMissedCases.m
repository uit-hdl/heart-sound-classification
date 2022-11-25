% plot AVPGmean, along with detected cases and missed cases. Run
% "jointVHDscreeningUsingAllPos" first. 

close
figure
plot(JvalTot,HSdata.avmeanpg(JvalTot),'k.')
hold on
plot(S.AS.J0.truePos,HSdata.avmeanpg(S.AS.J0.truePos),'.','color','g')
% plot(S.AS.J0.falsePos,HSdata.avmeanpg(S.AS.J0.falsePos),'d','color',color2triplet('red'),'MarkerSize',8)
plot(S.AS.J0.falsePos,HSdata.avmeanpg(S.AS.J0.falsePos),'r.')
plot(S.AS.J0.falseNeg,HSdata.avmeanpg(S.AS.J0.falseNeg),'.','color','red')
yline(15,'--')
xlim([0,max(JvalTot)+1])

