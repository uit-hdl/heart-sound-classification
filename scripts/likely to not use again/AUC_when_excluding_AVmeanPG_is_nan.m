% run script jointVHDscreeningUsingRiskFactorModel.m to obtain variables

I = isnan(HSdata.avmeanpg);
Npred = height(activations);
J = find( findInd(JvalTot, find(I)) );
J = setdiff(1:Npred, J);

activations(J)

[AUC,X,Y] = performanceSummaryNeurNet([],YtargetTot(J),activations(J),...
                                      [],[],plotValTot);

