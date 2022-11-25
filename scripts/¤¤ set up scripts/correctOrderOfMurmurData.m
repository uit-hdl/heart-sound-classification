% the numbering of the auscultation areas was changed during the collection
% of the data, and here I correct for that. The correction is defined by
% the map (1,2,3,4) -> (4,3,2,1). This mapping is performed on all rows
% that corresponds to dates that are earlier than Aug. 24 2015, as these
% follow the order 1:Mitral 2:tricuspid 3:Pulmonic 4:Aortic. I want to
% reverse the order of these.
HSdataCopy = HSdata;

echoDates = datetime(HSdata.ECHO_CHKIN_DATE_T72);
orderReversalDate = datetime(2015,08,24);

JbeforeW35 = find(echoDates < orderReversalDate);
JafterW35  = find(echoDates >= orderReversalDate);
idBeforeW35 = HSdata.UNIKT_LOPENR(JbeforeW35);
idAfterW35 = HSdata.UNIKT_LOPENR(JafterW35);

idWrong = idBeforeW35;
idCorr  = idAfterW35;
%%

Icorr  = findInd(HSdata.UNIKT_LOPENR,idCorr);
Iwrong = findInd(HSdata.UNIKT_LOPENR,idWrong);

%% find the groups of columns that you want to reverse the order for:

% observer AD, indicator Normal, systolic, diastolic, noise
ADind(1,:) = [41,45,49,53];
ADind(2,:) = ADind(1,:) + 1;
ADind(3,:) = ADind(2,:) + 1;
ADind(4,:) = ADind(3,:) + 1;

for i=1:4
    HSdataCopy(Iwrong,ADind(i,:)) = HSdataCopy(Iwrong,flip(ADind(i,:)));
end

% observer SA, indicator Normal, systolic, diastolic, noise
SAind(1,:) = [73,77,81,85];
SAind(2,:) = SAind(1,:) + 1;
SAind(3,:) = SAind(2,:) + 1;
SAind(4,:) = SAind(3,:) + 1;

for i=1:4
    HSdataCopy(Iwrong,SAind(i,:)) = HSdataCopy(Iwrong,flip(SAind(i,:)));
end

% murmur grade, observer AD
ADgrade = [57 58 59 60];
HSdataCopy(Iwrong,ADgrade) = HSdataCopy(Iwrong,flip(ADgrade));

% murmur grade, observer SA
SAgrade = [89 90 91 92];
HSdataCopy(Iwrong,SAgrade) = HSdataCopy(Iwrong,flip(SAgrade));

% audibility of second tone observer AD
% observer AD, indicator Normal, systolic, diastolic, noise
AD2tone(1,:) = 61+(0:3)*3;
AD2tone(2,:) = AD2tone(1,:) + 1;
AD2tone(3,:) = AD2tone(2,:) + 1;
for i=1:3
    cols = AD2tone(i,:);
    HSdataCopy(Iwrong,cols) = HSdataCopy(Iwrong,flip(cols));
end

% THERE IS A COLUMN MISSING HERE; SO ILL JUST LEAVE THESE ALONE
% audibility of second tone observer SA
% observer AD, indicator Normal, systolic, diastolic, noise
% SA2tone(1,:) = 93+(0:3)*3;
% SA2tone(2,:) = SA2tone(1,:) + 1;
% SA2tone(3,:) = SA2tone(2,:) + 1;
% for i=1:3
%     cols = SA2tone(i,:);
%     HSdataCopy(Iwrong,cols) = HSdataCopy(Iwrong,flip(cols));
% end


% sound agreement (binary) for each auscultation location
clear c
c(1,:) = [104 108 112 116];
c(2,:) = c(1,:) + 1;
c(3,:) = c(2,:) + 1;
c(4,:) = c(3,:) + 1;
for i=1:4
    cols = c(i,:);
    HSdataCopy(Iwrong,cols) = HSdataCopy(Iwrong,flip(cols));
end


% Faintness of first heart sound for each location, obs AD
clear c
c(1,:) = [122 124 126 128];
c(2,:) = c(1,:) + 1;
for i=1:2
    cols = c(i,:);
    HSdataCopy(Iwrong,cols) = HSdataCopy(Iwrong,flip(cols));
end

% Faintness of first heart sound for each location, obs SA
clear c
c(1,:) = [130 132 134 136];
c(2,:) = c(1,:) + 1;
for i=1:2
    cols = c(i,:);
    HSdataCopy(Iwrong,cols) = HSdataCopy(Iwrong,flip(cols));
end

%% All done; save modified data frame

HSdata = HSdataCopy;

