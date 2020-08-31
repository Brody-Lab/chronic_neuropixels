P=get_parameters;
if ~exist('Cells', 'var')
    load([P.data_folder_path filesep 'Cells.mat'])
    fprintf('\nLoaded the variabe CELLS.')
end
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
                                 'legend_on', k ~= 4, ...
                                 'ylabel_on', mod(k,5)==1)
    title(lookup_name_of_brain_area(brain_area), 'FontWeight', 'Normal')
    set(gca, 'Xtick', [1,4,16,64,256], ...
             'FontSize', 18)
    if k < 11, xlabel(''); end
end
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep mfilename], P.figure_image_format{i})
end