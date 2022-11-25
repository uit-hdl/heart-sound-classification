function I = id2ind(id,dataFrame)
% inverse function of ind2id which converts participant-id to its index I
% in dataFrame.
I = findInd(dataFrame.UNIKT_LOPENR,id);
end