P=get_parameters;
if ~exist('Cells', 'var')
    fprintf('\nLoading the variabe CELLS...')
    load([P.data_folder_path filesep 'Cells.mat'])
    fprintf('\nDone loading the variabe CELLS.')
end
%% Create a table
% the animal, age at time of implant, probe newness, probe tip depth,
% estimated number of electrodes recorded from in that region, and shank
% orientation; the (AP,ML,DV) coordinates of the approximate midpoint of the
% probe's location within that brain region

T = tabulate_implants_for_SU_curves(Cells);
writetable(T, [P.plots_folder_path, filesep 'figure2_supp_SU_curve.csv'])
%%
nrow = 4;
ncol = 5;
figure('Pos', [1,1,2560, 1440])

k = 0;
for brain_area = P.brain_areas(:)'
    k = k + 1;
    ax = subplot(nrow,ncol,k);
    plot_longevity_by_brain_area(Cells, char(brain_area), ...
                                 'ax', ax, ...
                                 'metric', 'single_unit', ...
                                 'normalize', true, ...
                                 'ylabel_on', k == 1, ...
                                 'T_implants', T)
    title(lookup_name_of_brain_area(brain_area))
    set(gca, 'Xtick', [1,4,16,64,256], ...
             'FontSize', 18)
    if k < 11, xlabel(''); end
end
saveas(gcf, [P.plots_folder_path, filesep 'figure2_supp_SU_curves'], 'svg')
