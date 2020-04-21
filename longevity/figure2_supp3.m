% Figure 2--supplement 3 fits a sum of two exponential decays to the unit
% counts
%
% The model:
%
%   n_units ~ alpha*exp(-(x-x0)/tau1) + (1-alpha)*exp(-(x-x0)/tau2) + epsilon_i
%
% The parameter "alpha" is the fraction of units that shows more rapid loss
% and is between 0 and 1. The parameter "tau1" is the time constant
% associated with the more rapid degradation, and "tau2" is associated with
% the slower time constant. The parameter "x0" is the earliest day since
% the surgery when we include data for analysis, and is 1. Finally,
% "epsilon_i" is i.i.d. Gaussians with the same variance. I have not
% determined whether a heteroskedastic model might change the conclusions.
%
% The unit counts for each dontion are normalized by the average on the
% first day after the surgery. This normalization removes one fewer
% parameter to estimate. 

% please specify me!
re_bootstrap = false;
%% Get files
if ~isfile(P.figure2_supp3_data_path)
    re_bootstrap = true;
end
P = get_parameters;
if re_bootstrap
    metric = 'unit';
    DV = struct;
    DV.T = get_metrics_from_Cells(Cells, 'condition_on', 'DV');
    [DV.p_hat, DV.p_CI, DV.p_boot] = fit_exp_decay_to_data(DV.T, ...
                                        'metric', metric, ...
                                        'normalize_initial_value', true);
    AP = struct;
    AP.T = get_metrics_from_Cells(Cells, 'condition_on', 'AP');
    [AP.p_hat, AP.p_CI, AP.p_boot] = fit_exp_decay_to_data(AP.T, ...
                                        'metric', metric, ...
                                        'normalize_initial_value', true);
    save(P.figure2_supp3_data_path, 'AP', 'DV', 'metric')
else
    load(P.figure2_supp3_data_path)
end
%%
figure('Pos', [100, 50, 1850, 700])
k = 0;
n_col = 4;
n_row = 2;
label_offset = 0;
ax_hdl = {};
plot_data = [true, false, false];
component = {'both', 'fast', 'slow'};
%% Different DV bins
for i = 1:3
    k = k + 1;
    ax_hdl{k} = subplot(n_row,n_col,k);
    plot_exp_decay_norm(DV.T, 'axes', ax_hdl{k}, ...
                                     'component', component{i}, ...
                                     'legend_on', i ==1, ...
                                     'metric', metric, ...
                                     'plot_data', plot_data(i), ...
                                     'xlabel', false, ...
                                     'ylabel', i ==1, ...
                                     'p_boot', DV.p_boot)
    label_panel(ax_hdl{k}, P.panel_labels(k+label_offset), ...
                'FontSize', P.panel_label_font_size);
end
% plot bootstrap samples
k = k + 1;
ax_hdl{k} = subplot(n_row,n_col,k);
set(gca, P.axes_properties{:})
set(gca, 'yscale', 'log')
for i=1:max(DV.T.condition)
    plot(DV.p_boot{i}(:,1), DV.p_boot{i}(:,3), '+', 'markersize', 3, 'Color', P.color_order(i,:))
end
for i=1:max(DV.T.condition)
    plot(DV.p_hat{i}(:,1), DV.p_hat{i}(:,3), ...
         'o', 'markersize', 10, ...
         'linewidth', 2, ...
         'Color', 'k')
end
ylabel('\tau_{slow}')
axes_pos_scaling = [1,1,1,1];
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0.01,0,0,0])
label_panel(ax_hdl{k}, P.panel_labels(k+label_offset), ...
                'FontSize', P.panel_label_font_size);
%% Different AP bins
clr_offset = max(DV.T.condition);
for i = 1:3
    k = k + 1;
    ax_hdl{k} = subplot(n_row,n_col,k);
    plot_exp_decay_norm(AP.T, 'axes', ax_hdl{k}, ...
                                     'component', component{i}, ...
                                     'legend_on', i ==1, ...
                                     'metric', metric, ...
                                     'plot_data', plot_data(i), ...
                                     'color_order_offset', clr_offset, ...
                                     'show_equation', false, ...
                                     'ylabel', i ==1, ...
                                     'p_boot', AP.p_boot)
    label_panel(ax_hdl{k}, P.panel_labels(k+label_offset), ...
                'FontSize', P.panel_label_font_size);
end

% plot bootstrap samples
k = k + 1;
ax_hdl{k} = subplot(n_row,n_col,k);
set(gca, P.axes_properties{:})
for i=1:max(AP.T.condition)
    plot(AP.p_boot{i}(:,1), AP.p_boot{i}(:,3), '+', 'markersize', 3, ...
        'Color', P.color_order(i+clr_offset,:))
    plot(AP.p_hat{i}(:,1), AP.p_hat{i}(:,3), ...
         'o', 'markersize', 10, ...
         'linewidth', 2, ...
         'Color', 'k')
end
set(gca, 'yscale', 'log')
h=xlabel('Fraction of units with fast decay (\alpha)', 'fontweight', 'normal')
ylabel('\tau_{slow}')
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0.01,0,0,0])
label_panel(ax_hdl{k}, P.panel_labels(k+label_offset), ...
                'FontSize', P.panel_label_font_size);
    %% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure2_supp3'], P.figure_image_format{i})
end