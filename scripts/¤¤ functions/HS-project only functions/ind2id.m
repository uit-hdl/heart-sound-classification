function id = ind2id(Index,dataFrame)
% Converts logical index vector or linear index which refers to rows in
% dataFrame, and returns a vector of the id's of those elements. The id's
% are given by the variable named UNIKT_LOPENR.
id = dataFrame.UNIKT_LOPENR(Index);

end