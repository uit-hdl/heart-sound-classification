%% 3d SCATTER PLOT
clear I
I{1} = HSdata.ARgrade==0;
I{2} = HSdata.ARgrade==1;
I{3} = HSdata.ARgrade==2;
I{4} = HSdata.ARgrade==3;
I{5} = HSdata.ARgrade==4;

xdata = HSdata.avarea;
ydata = HSdata.LVEFBIPLANE_T72;
zdata = HSdata.ARMAXPG_T72;

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
X{5} = xdata(I{5});
Y{5} = ydata(I{5});
Z{5} = zdata(I{5});

clf
scatter3(X{1},Y{1},Z{1},[],'g')
hold on
scatter3(X{2},Y{2},Z{2}, [],'b')
scatter3(X{3},Y{3},Z{3},[],'b')
scatter3(X{4},Y{4},Z{4},[],'r')
scatter3(X{5},Y{5},Z{5},[],'r')

% scatter3([1.5,1],[0,0],[0,0],[],'k','filled')
% scatter3([0,0],[3,4],[0,0],[],'k','filled')
% scatter3([0,0],[0,0],[25,40],[],'k','filled')

xlabel '1'
ylabel '2'
zlabel '3'