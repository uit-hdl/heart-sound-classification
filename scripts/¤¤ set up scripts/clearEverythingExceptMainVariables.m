% clear workspace variables that are probably not going to be used again:
initialVars = {'HMMpar' 'HSdata0' 'HSdata' ...
    'HSdataTrain' 'Jtest0' 'Jtrain0' 'Jval0' 'stats' 'statsAll' 'initialVars'};
clearvars('-except',initialVars{:})