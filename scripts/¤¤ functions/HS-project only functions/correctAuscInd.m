function y = correctAuscInd(id,auscIndex,idAuscOrder)
% this function takes the id of the subject and the auscultation index and
% maps it to the correct index (as in: 1-aortic, 2-pulmonic, 3-tricuspid,
% 4-mitral) by checking if the id is in the vector of id's that have
% reverse order. This needs to be done because after week 34, the order was
% reversed.

if ismember(id,idAuscOrder.APTM)
    y = auscIndex;
else
    y = 5 - auscIndex;
end
end