AP_bin_edges = [-8, 0, 4];
T = get_metrics_from_Cells(Cells, 'condition_on', 'AP', ...
                              'AP_bin_edges', AP_bin_edges);
%%
DV_bin_edges = [-10, -2,  -1,  0];
T = get_metrics_from_Cells(Cells, 'condition_on', 'DV', ...
                                  'DV_bin_edges', DV_bin_edges);
%%
[p_hat, p_CI, p_boot, exit_flag] = fit_exp_decay_to_data(T, 'normalize_initial_value', true, ...
                                                            
                                                 'n_boot', 1000)
%%
n_cond = max(T.condition);
P = get_parameters;

big_figure
subplot(3,3,1)
hold on
for i = 1:n_cond
    plot(p_boot{i}(:,1), p_boot{i}(:,2), 'o', 'Color', P.color_order(i,:))
end
set(gca,'yscale', 'log')
xlabel('frac')
ylabel('tau1')

subplot(3,3,2)
hold on
for i = 1:n_cond
    plot(p_boot{i}(:,1), p_boot{i}(:,3), 'o', 'Color', P.color_order(i,:))
end
set(gca,'yscale', 'log')
xlabel('frac')
ylabel('tau2')

subplot(3,3,3)
hold on
for i = 1:n_cond
    plot(p_boot{i}(:,2), p_boot{i}(:,3), 'o', 'Color', P.color_order(i,:))
end
set(gca, 'xscale', 'log', 'yscale', 'log')
xlabel('tau1')
ylabel('tau2')

subplot(3,3,4)
hold on
for i = 1:n_cond
    histogram(p_boot{i}(:,1), 0:0.04:1, 'FaceColor', P.color_order(i,:))
end
xlabel('frac')

subplot(3,3,5)
hold on
for i = 1:n_cond
    histogram(log(p_boot{i}(:,2)), 25, 'FaceColor', P.color_order(i,:))
end
xlabel('log(tau1)')

subplot(3,3,6)
hold on
for i = 1:n_cond
    histogram(log(p_boot{i}(:,3)), 25, 'FaceColor', P.color_order(i,:))
end
xlabel('log(tau2)')

