function [AUC,X,Y,T,p,YvalWhole,P] = performanceSummaryNeurNet(XvalSeg,YvalSeg,...
                                            net,Npcg,NsegPerPCG,figures,Nloc)
% Get AUC and plot sensitivity and specificity curve. Gets performance
% measures for both segmented data (treated as independent samples) and
% whole data, where each set of segment predictions gets combined into a
% single prediction using the median classification probability for each
% set. net is either the network or the activation output of the network on
% XvalSeg.

% Npcg = Number of in sample before
%% preliminary:
% Y is the sensitivity, and 1-X is the specificity
if nargin==5
    figures = true;
    Nloc = 1;
elseif nargin==6
    Nloc = 1;
end

if isempty(Npcg)
    % assume that the vectors provided corresponds to whole recordings
    Npcg = numel(YvalSeg);
end
if isempty(NsegPerPCG)
    NsegPerPCG = 6;
end
%% code:

if ~isnumeric(net)
    act = activations(net,XvalSeg,'softmax');
    p.mur = act(2,:);
elseif isnumeric(net)
    % if not categorical, then it is numerical
    if isvector(net)
        p.mur = net;
    else
        % the argument "net" is assumed to be a n by 2 matrix with
        % activations for each of the two classes:
        p.mur = net(:,2);
    end
    if iscategorical(YvalSeg)
        % convert to logical
        YvalSeg = YvalSeg=='1';
    end
end
% reshape back into form whith Nval rows and NsegPerPCG columns, to get label
% for each set of segments:
YvalWhole = mode(reshapeSegmentScores(YvalSeg,Npcg,NsegPerPCG,Nloc),2);

% p.murWhole = reshape(p.mur,[Nval,NsegPerPCG]);
p.murWhole = reshapeSegmentScores(p.mur,Npcg,NsegPerPCG,Nloc);
p.murWhole = median(p.murWhole,2);

if numel(YvalSeg)==Npcg && numel(p.mur)>Npcg
    % if validation labels contains signal labels and activations contain
    % segment predictions:
    YvalSeg = repmat(YvalSeg, 1,NsegPerPCG);
    YvalSeg = reshape(YvalSeg,[NsegPerPCG*Npcg,1]);
end

if sum(YvalSeg)>0
    [X.seg,Y.seg,T.seg,AUC.seg]         = perfcurve(YvalSeg, p.mur, true);
end
if sum(YvalWhole)>0
    [X.whole,Y.whole,T.whole,AUC.whole] = perfcurve(YvalWhole, p.murWhole, true);
end

if sum(YvalSeg)==0 || sum(YvalWhole)==0
    X = nan;
    Y = nan;
    T = nan;
    AUC.seg = nan;
    AUC.whole = nan;
end

% plot
if figures && sum(YvalWhole)>0
    % Y is the sensitivity, and 1-X is the specificity
    
    % ¤¤ CHOOSE DISPLAY STYLE ¤¤
    displayStyle = 2;
    if displayStyle==1
        plot(X.seg,Y.seg)
        hold on
        plot(X.whole, Y.whole)
        plot([0,0],[1,1],'g')
        title(sprintf('AUCseg=%.3g, AUCwhole=%.3g',AUC.seg,AUC.whole))
        legend({'segments','whole PCG-signal'})
    elseif displayStyle==2
        P = plot(X.whole,Y.whole);
%         hold on
%         plot([0,0],[1,1],'g')
%         title(sprintf('AUC=%.3g',AUC.whole))
    elseif displayStyle==3
        plot(X.whole, Y.whole)
        hold on
        plot([0,0],[1,1],'g')
        title ROC-curve
    end
    xlabel '1 - specificity'
    ylabel 'sensitivity'
end

end