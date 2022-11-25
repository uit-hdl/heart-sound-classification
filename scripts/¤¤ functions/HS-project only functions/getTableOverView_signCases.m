function T = getTableOverView_signCases(sigCaseMat,predMat,JvalTot)
% generates an overview (in the form of a matrix) of COMBINED SCREENING for
% SIGNIFICANT cases. sigCaseMat is a 4-column binary matrix (full form,
% with columns corresponds to AR,MR,AS, and MS, in that order) indicating
% the positions of the significant cases (either determined by grade, or
% grade and symptoms). Similarly, predMatrix is a binary 4 column matrix
% indicating predictions of significant cases for each VHD. JvalTot is a
% linear index vector that indicates in which rows are to be included in
% the analysis, and can be either prediction-rows, or rows corresponding to
% prediction and complete data.

% sig. case screening:
predAtleastOne = sum(predMat,2)>0;
YsigSickAtleastOne = sum(sigCaseMat,2)>0;

Ypred   = predAtleastOne(JvalTot);
Ytarget = YsigSickAtleastOne(JvalTot);

calcCI = true;

% sensitivity and specificity for detecting sig. cases:
[S.any.sn, sn_ci] =  condProb(Ypred,Ytarget,calcCI);
[S.any.sp, sp_ci] =  condProb(~Ypred,~Ytarget,calcCI);
S.any.sn_stdErr = (sn_ci(2)-sn_ci(1))/2;
S.any.sp_stdErr = (sp_ci(2)-sp_ci(1))/2;
S.any.Ncaught = sum(and(Ypred,Ytarget));
S.any.Nmissed = sum(and(~Ypred,Ytarget));
S.any.sn = presentFraction(S.any.sn);
S.any.sp = presentFraction(S.any.sp);
S.any.sn_stdErr = presentFraction(S.any.sn_stdErr);
S.any.sp_stdErr = presentFraction(S.any.sp_stdErr);

VHDnames = {'AR','MR','AS','MS'};
for i=1:4
    % SN and SP for specific diseases for sig. case screening:
    name = VHDnames{i};
    
    [S.(name).sn, sn_ci] =  condProb(Ypred,sigCaseMat(JvalTot,i),calcCI);
    [S.(name).sp, sp_ci] =  condProb(~Ypred,~sigCaseMat(JvalTot,i),calcCI);
    S.(name).sn_stdErr = (sn_ci(2)-sn_ci(1))/2;
    S.(name).sp_stdErr = (sp_ci(2)-sp_ci(1))/2;
    
    S.(name).Ncaught = sum(and(Ypred,sigCaseMat(JvalTot,i)));
    S.(name).Nmissed = sum(and(~Ypred,sigCaseMat(JvalTot,i)));
    
    S.(name).sn = presentFraction(S.(name).sn);
    S.(name).sp = presentFraction(S.(name).sp);
    S.(name).sn_stdErr = presentFraction(S.(name).sn_stdErr);
    S.(name).sp_stdErr = presentFraction(S.(name).sp_stdErr);
end

% summarize info in table:
T = [struct2table(S.AR);...
     struct2table(S.MR);...
     struct2table(S.AS);...
     struct2table(S.MS);...
     struct2table(S.any)];

T.Properties.RowNames = {'AR','MR','AS','MS','allVHD'};
T = rows2vars(T);

end