% I run this after setting the date column to increasing order (dates are increasing)
idAuscOrder.MTPA = HSdata0.UNIKT_LOPENR(1:29);
idAuscOrder.APTM = HSdata0.UNIKT_LOPENR(30:height(HSdata0));
save('idAuscOrder.mat','idAuscOrder')