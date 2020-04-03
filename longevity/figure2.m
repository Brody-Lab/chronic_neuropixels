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
figure('Pos', [100, 50, 2000, 1300])
k = 0;
n_col = 4;
n_row = 4;
label_offset = 0;
% -------------------------------------------------------------------------
% Plot individual recordings
for region = {'mFC', 'MCtx_ADS', 'AVS'}
    k = k + 1;
    ax_hdl(k) = subplot(n_row, n_col,k);
    plot_indiv_recordings(Cells, region{:}, ...
                          'metric', 'n_good_units', ...
                          'axes', gca,...
                          'ylabel_on', mod(k, n_col)==1);
end
% shift the plots right and above
inter_axes_dx = ax_hdl(1).Position(1) + ax_hdl(2).Position(1);
dx = inter_axes_dx/(n_col-1)/2;
for i = 1:3
    ax_hdl(i).Position(1)=ax_hdl(i).Position(1)+dx;
    ax_hdl(i).OuterPosition(2) = ax_hdl(i).OuterPosition(2) + 0.03;
    label_panel(ax_hdl(i), P.panel_labels(i), 'FontSize', P.panel_label_font_size);
end
% -------------------------------------------------------------------------
% move on to the next row
k=k+1;
label_offset = -1;
% -------------------------------------------------------------------------
% Plot average across all conditions
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k +     1;
    subplot(n_col, n_row,k);
    plot_average_longevity(Cells, 'metric', metric{:}, ...
                                  'axes', gca);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    if mod(k,n_col)==1
        legend('Overall average')
    end
    title(P.text.(metric{:}))
end
% -------------------------------------------------------------------------
% Plot average, conditioned on DV
bin_edges = [-10, -2,  -1,  0];
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_col, n_row,k);
    legend_on = mod(k,n_col)==1;
    plot_average_longevity(Cells, 'metric', metric{:}, ...
                                  'condition_on', 'DV', ...
                                  'anatom_bin_edges', bin_edges, ...
                                  'legend_on', legend_on, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3, ...
                                  'axes', gca);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    
    % Report the stats from the ANOVA
    Pval.(metric{:}) = examine_position_dependence(Cells, 'condition_on', 'DV', ...
                                                     'anatom_bin_edges', bin_edges, ...
                                                     'elect_bin_edges', [1, 385, 769], ...
                                                     'metric', metric{:});
    fprintf('\nANOVA with factors 1)bank and 2)DV on the change in %s:', metric{:});
    fprintf('\n    bank: p = %0.3f', Pval.(metric{:})(1));
    fprintf('\n    DV: p = %0.3f', Pval.(metric{:})(2));
end
% -------------------------------------------------------------------------
% Plot average, conditioned on AP
bin_edges = [-8, 0, 4];
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}
    k = k + 1;
    subplot(n_col, n_row,k);
    legend_on = mod(k,n_row)==1;
    plot_average_longevity(Cells, 'metric', metric{:}, ...
                                  'color_order_offset', 3, ...
                                  'condition_on', 'AP', ...
                                  'anatom_bin_edges', bin_edges, ...
                                  'legend_on', legend_on, ...
                                  'normalize_by_electrodes', true, ...
                                  'FaceAlpha', 0.3, ...
                                  'axes', gca);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
end
% -------------------------------------------------------------
% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure2'], P.figure_image_format{i})
end
