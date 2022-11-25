clear I
I{1} = HSdata.ASgrade==0;
I{2} = HSdata.ASgrade==1;
I{3} = HSdata.ASgrade==2;
I{4} = HSdata.ASgrade==3;

xdata = HSdata.avmeanpg;
ydata = HSdata.AVMAXPG_T72;

X{1} = xdata(I{1});
Y{1} = ydata(I{1});
X{2} = xdata(I{2});
Y{2} = ydata(I{2});
X{3} = xdata(I{3});
Y{3} = ydata(I{3});
X{4} = xdata(I{4});
Y{4} = ydata(I{4});

clf
scatter(X{1},Y{1},[],'g')
hold on
scatter(X{2},Y{2},[],'b')
scatter(X{3},Y{3},[],'k')
scatter(X{4},Y{4},[],'r')
xlabel 'cardiac output'
ylabel 'mean pressure gradient' 

%% 3d SCATTER PLOT

I{1} = stats.I.AS.EQ.grade0.any;
I{2} = stats.I.AS.EQ.grade1.any;
I{3} = stats.I.AS.EQ.grade2.any;
I{4} = stats.I.AS.EQ.grade3.any;

xdata = HSdataTrainAndVal.AVAVMAX_T72;
ydata = log(HSdataTrainAndVal.AVVMAX_T72);
zdata = HSdataTrainAndVal.AVMEANPG_T72;

X{1} = xdata(I{1});
Y{1} = ydata(I{1});
Z{1} = zdata(I{1});
X{2} = xdata(I{2});
Y{2} = ydata(I{2});
Z{2} = zdata(I{2});
X{3} = xdata(I{3});
Y{3} = ydata(I{3});
Z{3} = zdata(I{3});
X{4} = xdata(I{4});
Y{4} = ydata(I{4});
Z{4} = zdata(I{4});

clf
scatter3(X{1},Y{1},Z{1},[],'g')
hold on
scatter3(X{2},Y{2},Z{2}, [],'b')

scatter3(X{3},Y{3},Z{3},[],'k')
scatter3(X{4},Y{4},Z{4},[],'r')
scatter3([1.5,1],[0,0],[0,0],[],'k','filled')
scatter3([0,0],[3,4],[0,0],[],'k','filled')
scatter3([0,0],[0,0],[25,40],[],'k','filled')

xlabel 'AV-area'
ylabel 'AV-velocity'
zlabel 'AV-pressure gradient'
