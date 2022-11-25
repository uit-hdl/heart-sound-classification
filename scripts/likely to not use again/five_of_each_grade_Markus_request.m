clear J
J{1} = find(HSdata.murGrade1==1,5);
J{2} = find(HSdata.murGrade1==2,5);
J{3} = find(HSdata.murGrade1==3,5);
J{4} = find(HSdata.murGrade1==4,5);
J{5} = find(HSdata.murGrade1==5,5);

clear AVmeanPG
fileNames = cell(5,5);
for grade=1:5
    for j=1:5
        id = HSdata.id(J{grade}(j));
        fileNames{grade,j} = sprintf('%.0f_hjertelyd_%g.wav',id,1);
        s = sprintf('grade%g',grade)
        AVmeanPG.(s).avmeanpg(j,:) = HSdata.avmeanpg(HSdata.id==id);
        AVmeanPG.(s).id(j,:) = string(id);
    end
end


%%
[x,fs] = wav2TS(10257520,1);
Nds = 20;
x = downsample(x,Nds);
fs = fs/Nds;
scaleogramPlot(x,fs);
%%
load networksTrainingSet_valStop.mat
predictMurmur(id,net)

