% FIGURE2 Make figure 2 from the manuscrip written by Luo, Bondy, et al.
%
%   The figure show that after an initial loss of units, spiking signals
%   can be maintained >60 days in anterior, deeper brain regions.

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
figure('Pos', [100, 50, 2000, 1300])
k = 0;
n_col = 4;
n_row = 4;
label_offset = 0;
show_fit = true;
ax_hdl={};
%% Plot individual recordings
for region = {'mFC', 'MCtx_ADS', 'AVS'}
    k = k + 1;
    ax_hdl{k} = subplot(n_row, n_col,k);
    plot_indiv_recordings(Cells, region{:}, ...
                          'metric', 'n_good_units', ...
                          'axes', gca,...
                          'ylabel_on', mod(k, n_col)==1);
end
% shift the plots right and above
inter_axes_dx = ax_hdl{1}.Position(1) + ax_hdl{2}.Position(1);
dx = inter_axes_dx/(n_col-1)/2;
for i = 1:3
    ax_hdl{i}.Position(1)=ax_hdl{i}.Position(1)+dx;
    ax_hdl{i}.OuterPosition(2) = ax_hdl{i}.OuterPosition(2) + 0.03;
    label_panel(ax_hdl{i}, P.panel_labels(i), 'FontSize', P.panel_label_font_size);
end
% move on to the next row
k=k+1;
label_offset = -1;
%% Plot average across all conditions
T = get_metrics_from_Cells(Cells);
fprintf('\nANOVA with the factor days since implants and for data > 7 days:')
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k +     1;
    subplot(n_col, n_row,k);
    if any(strcmp(metric{:}, {'unit', 'single_unit'})) && show_fit
        fit_type = 'exponential';
    else
        fit_type='';
    end
    plot_average_stability(T, 'metric', metric{:}, ...
                              'axes', gca, ...
                              'print_sample_size', mod(k,n_col)==1, ...
                              'fit_type', fit_type);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    title(P.text.(metric{:}))
    % report the pval based on the t-statistic
    pval = stat_test_stability(T, 'metric', metric{:}, 'day_range', [8, inf]);
    fprintf('\n    %s: p=%0.3f', metric{:}, pval);
end
fprintf('\n')
%% Plot average, conditioned on DV
T = get_metrics_from_Cells(Cells, 'condition_on', 'DV');
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_col, n_row,k);
    if any(strcmp(metric{:}, {'unit', 'single_unit'})) && show_fit
        fit_type = 'exponential';
    else
        fit_type='';
    end
    plot_average_stability(T, 'metric', metric{:}, ...
                              'fit_type', fit_type, ...
                                  'legend_on', false, ...
                                  'print_sample_size', mod(k,n_col)==1, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3, ...
                                  'axes', gca);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    % Report the stats from the ANOVA
    Pval.(metric{:}) = anova_DV_shank(Cells, 'DV_bin_edges', P.DV_bin_edges, ...
                                              'EI_bin_edges', P.EI_bin_edges, ...
                                                     'metric', metric{:});
    fprintf('\nANOVA with factors 1)bank and 2)DV on the change in %s:', metric{:});
    fprintf('\n    bank: p = %0.3f', Pval.(metric{:})(1));
    fprintf('\n    DV: p = %0.3f', Pval.(metric{:})(2));
end
%% Plot average, conditioned on AP
AP_bin_edges = [-8, 0, 4];
T = get_metrics_from_Cells(Cells, 'condition_on', 'AP', ...
                                  'AP_bin_edges', AP_bin_edges);
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_col, n_row,k);
    legend_on = mod(k,n_row)==1;
    if any(strcmp(metric{:}, {'unit', 'single_unit'})) && show_fit
        fit_type = 'exponential';
    else
        fit_type='';
    end
    plot_average_stability(T, 'metric', metric{:}, ...
                                 'fit_type', fit_type, ...
                                  'color_order_offset', 3, ...
                                  'legend_on', false, ...
                                  'print_sample_size', mod(k,n_col)==1, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3, ...
                                  'axes', gca);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
end
%% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure2_v1'], P.figure_image_format{i})
end