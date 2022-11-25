function [activ,u0,Ytarget,Ypred,AUC,glm] = get_sigVHDactivations(ActMat,dataFrame,Itrain,...
                                                    Ival,disease,classThr,plotVal,...
                                                    minSn,minSp,predWithMaxGrade)
% computes (joint) activation vector for the validation set. How the
% activation vector is computed depends on which disease is to be
% predicted. AS is predicted using the AVmeanPG calibration method, whereas
% the other diseases simply uses the maximum activation across the
% positions. ActMat is a 4-by-Ndata (2124 for HSdata) zero-padded matrix
% with murmur-algorithm activations, and Itrain and Ival are index vectors
% that gives index for the positions for which prediction is possible, i.e.
% positions where atleast one position is associated with a prediction.
if nargin==6
    plotVal = false;
    minSn = 0;
    minSp = 0;
    predWithMaxGrade = false;
elseif nargin==7
    minSn = 0;
    minSp = 0;
    predWithMaxGrade = false;
elseif nargin==9
    predWithMaxGrade = false;
end

if disease=="AS" && ~predWithMaxGrade
    [activ,u0,Ytarget,Ypred,AUC,glm] = get_ASactivations(ActMat,dataFrame,Itrain,Ival,...
                                                classThr,plotVal,minSn,minSp);
else
    A = max(ActMat,[],2);
    
    % *** estimate optimal threshold ***
    if length(disease)>2
        varName = disease;
    else
        varName = sprintf('%sgrade',disease);
    end
    Ytarget.train = dataFrame.(varName)(Itrain) >= classThr;
    activ.train   = A(Itrain);
    [~,X,Y,T,~] = performanceSummaryNeurNet([],Ytarget.train,activ.train,...
                                            [],[],false); 
    
    u0 = getOptimalThr(X,Y,T,minSn,minSp);
    
    % *** performance on validation set ***
    Ytarget.val = dataFrame.(varName)(Ival) >= classThr;
    activ.val   = A(Ival);
    AUC = performanceSummaryNeurNet([],Ytarget.val,activ.val,...
                                            [],[],false);   
    AUC = AUC.whole;
    
    Ypred.val = activ.val>=u0;
    
    glm = [];
end
end