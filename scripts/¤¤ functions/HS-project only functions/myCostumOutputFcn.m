function stop = myCostumOutputFcn(info)

stop = false;

persistent bestValAccuracy
persistent currentEpoch
persistent valLossHistory
persistent stoppageTicker


if info.State == "start"
    bestValAccuracy = 0;
    currentEpoch    = 1;
    valLossHistory = [];
    stoppageTicker = 0;
end

if ~isempty(info.ValidationRMSE)
    currentEpoch = info.Epoch;
    bestValAccuracy = max(bestValAccuracy,info.ValidationAccuracy);
    valLossHistory = [valLossHistory,info.ValidationLoss];

    if currentEpoch>5
        valLossHistory = valLossHistory(2:end);
        valLossMin = min(valLossHistory);
        noValImprovement = valLossMin<valLossHistory(end);
        if noValImprovement
            stoppageTicker = stoppageTicker + noValImprovement;
        else
            stoppageTicker = 0;
        end
    end        
    
    display(info.ValidationLoss)
    display(stoppageTicker)
    display(valLossHistory)
    
    if stoppageTicker>=7
        stop = true;
    end

end

end