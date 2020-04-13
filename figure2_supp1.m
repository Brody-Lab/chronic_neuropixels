% FIGURE2_SUPP1 Make figure 2--supplement 1 from the manuscrip written by
% Luo, Bondy, et al.
%
%   The figure show that no degradation of spiking signals was detected
%   over two months in rat medial prefrontal cortex (mPFC), the brain
%   region in which the stability of spiking signals was examined in (Jun
%   et al., 2017)
%% Load Cells
from_scratch = false; % Do you want to reassemble the data files from scratch?
P=get_parameters;
if from_scratch
    collect_cell_files
    postprocess_Cells
else
    if ~exist('Cells', 'var')
        fprintf('Loading the variabe CELLS...')
        load([P.data_folder_path filesep 'Cells.mat'])
    end
end
%% Make the figue
figure('Pos', [100, 50, 2000, 300])
k = 0;
n_col = 4;
n_row = 1;
label_offset = 0;
T = get_metrics_from_Cells(Cells, 'condition_on', 'brain_area', ...
                                  'brain_area', {{'PrL', 'MO'}, {'other'}});
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_row, n_col,k);
    plot_average_stability(T, 'metric', metric{:}, ...
                              'axes', gca, ...
                              'fit_type', '');
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    title(P.text.(metric{:}))
end
%% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure_compare_mPFC'], P.figure_image_format{i})
end