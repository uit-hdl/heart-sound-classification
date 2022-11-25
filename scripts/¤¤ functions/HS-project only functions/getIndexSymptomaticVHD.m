function sympSigVHD = getIndexSymptomaticVHD(dataFrame,classThr,symptomNames)
% gets index of clinically relevant disease, defined as symptomatic AR or
% MR grade 3 or higher, and presence of stenosis. symptomNames containes
% the names of the symptoms of AR and MR.

if nargin==2
    symptomCellArray = dataFrame.anginaOrDyspnea;
else

    Nsymp = length(symptomNames);
    symptomCellArray = cell(1,Nsymp);
    for i=1:Nsymp
        symptomCellArray{i} = dataFrame.(symptomNames{i})>0;
    end
end


IconfirmedARsymp = myor(symptomCellArray);
confirmedMRsymp = myor(symptomCellArray);
confirmedASsymp = dataFrame.ASgrade>=classThr(3);
confirmedMSsymp = dataFrame.MSgrade>=classThr(4);

sympSigVHD(:,1) = and(IconfirmedARsymp,dataFrame.ARgrade>=classThr(1));
sympSigVHD(:,2) = and(confirmedMRsymp,dataFrame.MRgrade>=classThr(2));
sympSigVHD(:,3) = confirmedASsymp;
sympSigVHD(:,4) = confirmedMSsymp;
sympSigVHD(:,5) = sum(sympSigVHD,2)>0;
     

end