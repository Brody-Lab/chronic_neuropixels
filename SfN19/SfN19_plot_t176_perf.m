fig_hdl = figure;

Pharma = PB_import_pharmacology_log;
idx = Pharma.Rat == 'T176' & (Pharma.Manipulation == 'saline' | (Pharma.Manipulation == 'muscimol' & Pharma.Hemisphere == 'bilateral'));
PB_plot_performance('T176', Pharma.Date(idx), 'ax_user', gca)
%%
title('')
legend off
set(gca, 'YGrid', 'off', ...
         'YTick', 0:20:100)
% text(-35, 70, 'saline', 'fontsize', 16)
% text(-35, 95, 'muscimol in SC', 'fontsize', 16)
set(gcf, 'Position', [ 1000          1094.33333333333                       261          243.66666666666])
saveas(gcf, 'C:\ratter\Analysis\tzluo\Plots\SfN19_plot_T176_perf\behav', 'svg')