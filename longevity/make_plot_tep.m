%figure2_v2_glm
tic
S = fit_exp_to_electrodes(Cells);
toc
b =[S.b];
%%
n_r = numel(P.exp_decay_regressors);
n = n_r+1;
P = get_parameters;
figure('Pos', [100, 50, 1500, 800])
k = 0;
n_row = 2;
n_col = 3;
label_offset = 8;
ax_hdl={};
XTickLabelRotation = 0;
k = 4;
axes_pos_scaling = [1,1,1,0.7];
axes_pos_offset = [0,0.04,0,0];

% Initial value
k = k + 1;
ax_hdl=subplot(n_row,n_col, k);
xticklabel = cellfun(@(x) P.text.(x), P.exp_decay_regressors, 'uni', 0);
beta_text = cellfun(@(x) ['e^{\beta' num2str(x) '}'], num2cell(1:n_r), 'uni', 0);
xticklabel = cellfun(@(x,y) [x,' \newline(' y ')'], xticklabel, beta_text, 'uni', 0);
set(gca, P.axes_properties{:}, ...
         'XLim', 0.5 + [1,n], ...
         'XTick', 2:n, ...
         'XTickLabel', xticklabel, ...
         'XTickLabelRotation', XTickLabelRotation)
plot(xlim, [0,0], 'k-')
b0 = b(2:n,:)';
med = median(b0);
CI = quantile(b0, [0.025, 0.975]);
lower = med - CI(1,:);
upper = CI(2,:) - med;
errorbar(2:n, med, lower, upper, 'ko',  'linewidth', 1)
title('Multipl. weight on the init. unit count', 'fontweight', 'normal')
ylabel('Weight',  'fontweight', 'normal');
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+axes_pos_offset)
label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);

% additive weight on the decay rate
k = k + 1;
subplot(n_row,n_col, k);
xticklabel = cellfun(@(x) P.text.(x), P.exp_decay_regressors, 'uni', 0);
beta_text = cellfun(@(x) ['\beta_{' num2str(x) '}'], num2cell((n_r+2):((n_r+1)*2-1)), 'uni', 0);
xticklabel = cellfun(@(x,y) [x,' \newline(' y ')'], xticklabel, beta_text, 'uni', 0);
ax_hdl=set(gca, P.axes_properties{:}, ...
         'XLim', 0.5 + [1,n], ...
         'XTick', 2:n, ...
         'XTickLabel', xticklabel, ...
         'XTickLabelRotation', XTickLabelRotation);
plot(xlim, [0,0], 'k-')
b1 = b(n_r+3:end,:)';
med = median(b1);
CI = quantile(b1, [0.025, 0.975]);
lower = med - CI(1,:);
upper = CI(2,:) - med;
errorbar(2:n, med, lower, upper, 'ko', 'linewidth', 1)
title('Additive weight on the decay rate', 'fontweight', 'normal')
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+axes_pos_offset);
label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);

% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure2_v2_glm'], P.figure_image_format{i})
end