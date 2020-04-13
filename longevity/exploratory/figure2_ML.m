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
figure('Pos', [100, 50, 2000, 300])
k = 0;
n_col = 4;
n_row = 1;
label_offset = 0;
% -------------------------------------------------------------------------
% Plot average across all conditions
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_row, n_col,k);
    plot_average_longevity(Cells, 'metric', metric{:}, ...
                                  'axes', gca, ...
                                  'normalize', true, ...
                                  'legend_on', mod(k,n_col)==1, ...
                                  'anatom_bin_edges', [0, 2, 4, 6], ...
                                  'condition_on', 'ML');
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    title(P.text.(metric{:}))
end
% Save
% for i = 1:numel(P.figure_image_format)
%     saveas(gcf, [P.plots_folder_path filesep 'figure_compare_mPFC'], P.figure_image_format{i})
% end