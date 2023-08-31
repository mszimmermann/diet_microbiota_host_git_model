function [gof] = calculate_gof(y_predicted, y_actual)

%normalize predicted rdata to max
cur_rdata = (y_predicted-0.8);
cur_rdata = cur_rdata/max(max(cur_rdata));

y_predicted = cur_rdata;

SSR = sum((y_predicted - y_actual).^2);
TSS = sum((y_actual - mean(y_actual)).^2);

Rsquared = 1 - SSR/TSS;

gof = Rsquared;