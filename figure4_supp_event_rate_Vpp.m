%% Load the cells data
from_scratch = false; % Do you want to reassemble the data files from scratch?
P=get_parameters;
if from_scratch
    Cells=collect_cells_files();
    Cells=postprocess_Cells(Cells);
else
    if ~exist('Cells', 'var')
        fprintf('Loading the variabe CELLS...')
        load([P.data_folder_path filesep 'Cells.mat'])
    end
end
%% Set up the figue
figure('Pos', [100, 50, 1455, 900])
k = 0;
n_row = 3;
n_col = 3;
label_offset = 0;
show_fit = true;
ax_hdl={};
%% A,B,C - Plot average across all conditions
fprintf('\nANOVA with the factor days since implants and for data > 7 days:')
for metric = {'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_row, n_col,k);
    exclude_multiunit = strcmp(metric, 'Vpp');
    T = get_metrics_from_Cells(Cells, 'exclude_multiunit', exclude_multiunit);
    plot_average_stability(T, 'metric', metric{:}, ...
                              'axes', gca, ...
                              'normalize_by_electrodes', true, ...
                              'print_sample_size', mod(k,n_row)==1);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    title(P.text.(metric{:}))
    % report the pval based on the t-statistic
    pval = stat_test_stability(T, 'metric', metric{:}, 'day_range', [8, inf]);
    fprintf('\n    %s: p=%0.3f', metric{:}, pval);
end
fprintf('\n')

% plot histogram
k = k + 1;
subplot(n_row, n_col,k)
T = get_metrics_from_Cells(Cells, 'exclude_multiunit', true);
make_Vpp_histogram(T)
label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size, 'x', -0.02);
%% D, E, F - Plot average, conditioned on DV
for metric = {'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_row, n_col,k);
    legend_on = mod(k,n_col)==1;
    exclude_multiunit = strcmp(metric, 'Vpp');
    T = get_metrics_from_Cells(Cells, 'condition_on', 'DV', ...
                                      'exclude_multiunit', exclude_multiunit);
    plot_average_stability(T, 'metric', metric{:}, ...
                                  'legend_on', legend_on, ...
                                  'print_sample_size', mod(k,n_row)==1, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3, ...
                                  'axes', gca);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
end

% plot histogram
k = k + 1;
subplot(n_row, n_col,k)
T = get_metrics_from_Cells(Cells, 'condition_on', 'DV', ...
                                  'exclude_multiunit', true);
make_Vpp_histogram(T)
label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size, 'x', -0.02);
%% G,H, I - Plot average, conditioned on AP
for metric = {'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_row, n_col,k);
    legend_on = mod(k,n_col)==1;
    exclude_multiunit = strcmp(metric, 'Vpp');
    T = get_metrics_from_Cells(Cells, 'condition_on', 'AP', ...
                                      'exclude_multiunit', exclude_multiunit);
    plot_average_stability(T, 'metric', metric{:}, ...
                                  'color_order_offset', 3, ...
                                  'legend_on', false, ...
                                   'legend_on', legend_on, ...
                                  'print_sample_size', mod(k,n_row)==1, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3, ...
                                  'axes', gca);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
end

% plot histogram
k = k + 1;
subplot(n_row, n_col,k)
T = get_metrics_from_Cells(Cells, 'condition_on', 'AP', ...
                                  'exclude_multiunit', true);
make_Vpp_histogram(T, 'color_order_offset', 3)
label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size, 'x', -0.02);
%% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep mfilename], P.figure_image_format{i})
end