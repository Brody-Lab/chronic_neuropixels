function[]=make_exp_decay_figure(T_ed, varargin)


%% Figure
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
% set(gca, 'yscale', 'log')
% for i=1:max(DV.T.condition)
%     plot(DV.p_boot{i}(:,1), DV.p_boot{i}(:,3), '+', 'markersize', 3, 'Color', P.color_order(i,:))
% end
% for i=1:max(DV.T.condition)
%     plot(DV.p_hat{i}(:,1), DV.p_hat{i}(:,3), ...
%          'o', 'markersize', 10, ...
%          'linewidth', 2, ...
%          'Color', 'k')
% end
% ylabel('\tau_{slow}')
% axes_pos_scaling = [1,1,1,1];
% set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0.01,0,0,0])
% label_panel(ax_hdl{k}, P.panel_labels(k+label_offset), ...
%                 'FontSize', P.panel_label_font_size);
syms x
for i=1:max(DV.T.condition)
    for j = 1:size(DV.p_boot{i},1)
        
        alpha = DV.p_boot{i}(j,1);
        tau_f = DV.p_boot{i}(j,2);
        tau_s = DV.p_boot{i}(j,3);
        
        eqn = exp(-1) == alpha*(exp(-1/tau_f))^(x-1) + (1-alpha)*(exp(-1/tau_s))^(x-1);
        
        tau_e{i}(j,1) = vpasolve(eqn, x);
        fprintf('%i',j)
    end
end
tau_e_med = cellfun(@median, tau_e);
tau_e_lower = tau_e_med - cellfun(@(x) quantile(x, 0.025), tau_e);
tau_e_upper = cellfun(@(x) quantile(x, 0.025), tau_e)-tau_e_med;

errorbar(1:max(DV.T.condition), tau_e_med, tau_e_lower, tau_e_upper);
set(gca, 'xLim', [0.5 max(DV.T.condition)+0.5], 'yscale', 'log')
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
% for i=1:max(AP.T.condition)
%     plot(AP.p_boot{i}(:,1), AP.p_boot{i}(:,3), '+', 'markersize', 3, ...
%         'Color', P.color_order(i+clr_offset,:))
%     plot(AP.p_hat{i}(:,1), AP.p_hat{i}(:,3), ...
%          'o', 'markersize', 10, ...
%          'linewidth', 2, ...
%          'Color', 'k')
% end
% set(gca, 'yscale', 'log')
% h=xlabel('Fraction of units with fast decay (\alpha)', 'fontweight', 'normal')
% ylabel('\tau_{slow}')
% set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0.01,0,0,0])
% label_panel(ax_hdl{k}, P.panel_labels(k+label_offset), ...
%                 'FontSize', P.panel_label_font_size);
tau_avg={};     
for i=1:max(AP.T.condition)
    tau_avg{i} = AP.p_boot{i}(:,1)./AP.p_boot{i}(:,2) + ...
                   (1-AP.p_boot{i}(:,1))./AP.p_boot{i}(:,3);
end
tau_e_med = cellfun(@median, tau_avg);
tau_e_lower = tau_e_med - cellfun(@(x) quantile(x, 0.025), tau_avg);
tau_e_upper = cellfun(@(x) quantile(x, 0.025), tau_avg)-tau_e_med;

errorbar(1:max(AP.T.condition), tau_e_med, tau_e_lower, tau_e_upper);
set(gca, 'xLim', [0.5 max(AP.T.condition)+0.5], 'yscale', 'log')

    %% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure2_supp3'], P.figure_image_format{i})
end