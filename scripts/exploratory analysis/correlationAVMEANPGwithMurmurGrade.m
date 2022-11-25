modelData = HSdata(union(Jtrain0,Jval0),:);
%#ok<*NOPTS>

InotNan = ~isnan(modelData.AVMEANPG_T72);
modelData = modelData(InotNan,:);

clear T
Tcorr = zeros(3,5);
Tauc  = zeros(5,5); 
for aa=1:5
    if aa<5
        murName = sprintf('Murmur_%g_grade_ref_ny_T72',aa);
    else
        murName = 'maxMeanMurGrade';
    end
    mycorr = @(x1,x2) corr(x1,x2);
    nIterations = 1000;
    % compute the CI of correlation between murmur grade and
    % aortic-valve-pressure-gradient:
    CI = bootci(nIterations,...
        {mycorr, modelData.(murName), modelData.AVMEANPG_T72} );

    Tcorr(:,aa) = [mycorr(modelData.(murName),modelData.AVMEANPG_T72);CI];
    for i=1:5
        [X,Y,~,A] = perfcurve(modelData.AVMEANPG_T72>=1+4*i,modelData.(murName),true); %#ok<*ASGLU>
        Tauc(i,aa) = A;
    end
end

modelData = HSdata;
InotNan = ~isnan(modelData.ASGRADE_T72);
modelData = HSdata(InotNan,:);

corrTable = array2table(Tcorr,'Var',{'pos1','pos2','pos3','pos4','max'},...
                'R',{'est.','ci lower','ci upper'})
AUCtable = array2table(Tauc,'Var',{'pos1','pos2','pos3','pos4','max'},...
                'R',{'thr. 5','thr. 9','thr. 14','thr. 18','thr. 22'})

[X,Y,~,A] = perfcurve(modelData.ASGRADE_T72>=1,modelData.AVMEANPG_T72,true);
% conclusion: AS and mean-pressure-gradient are near the identical.
[X,Y,~,A] = perfcurve(modelData.AVMEANPG_T72>=15,modelData.maxMeanMurGrade,true);
plot(X,Y)
title(sprintf('ROC curve MMG predicting AVMEANPG, AUC=%.3g',A*100))
% conclusion: maxMurmurGrade is higly predictive of mean-pressure-gradient.
%%
mycorr = @(x1,x2) corr(x1,x2);
nIterations = 1000;
CI = bootci(nIterations,...
    {mycorr,modelData.Murmur_1_grade_ref_ny_T72,modelData.ASGRADE_T72});
mycorr(modelData.Murmur_1_grade_ref_ny_T72,modelData.ASGRADE_T72)
CI'

