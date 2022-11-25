% What happens when we use cutoff murmur-grade 1,2 and 3? To get required
% variables, run: jointVHDscreeningUsingAllPos
%#ok<*NOPTS>
SN = zeros(3,3);
SP = zeros(3,3);
XX = cell(3,1);
YY = cell(3,1);
%%
murThr = [0.5 1 2];
for u=1:3
if targetType=="sigVHD31"
    target = 1;
elseif targetType=="sigVHD41"
    target = 2;
elseif targetType=="sigSymptVHD"
    target = 3;
end
    
Ypred = activations>murThr(u);
Ytarget = YtargetTot;

[~,X,Y] = performanceSummaryNeurNet([],YtargetTot,activations,...
                              [],[],plotValTot);
SN(target,u) = condProb(Ypred,Ytarget)
SP(target,u) = condProb(~Ypred,~Ytarget)
end
XX{target} = X.whole 
YY{target} = Y.whole
                          



%%
close all
plot(XX{1},YY{1})
hold on
plot(XX{2},YY{2})
plot(XX{3},YY{3})
xlabel '1-specificity'
ylabel 'sensitivity'

plot(1-SP(:,1),SN(:,1),'o','color','k','MarkerFaceColor','c')
plot(1-SP(:,2),SN(:,2),'o','color','k','MarkerFaceColor','g')
plot(1-SP(:,3),SN(:,3),'o','color','k','MarkerFaceColor','b')
legend('AR>=3, MR>=3, AS>=1 or MS>=1 (AUC=0.61)',...
    'AR>=4, MR>=4, AS>=1 or MS>=1 (AUC=0.731)',...
    'symptomatic AR>=3, MR>=3, or any stenosis (AUC=0.882)',...
    'murmur-grade-threshold = 0.5','murmur-grade-threshold = 1','murmur-grade-threshold = 2')






