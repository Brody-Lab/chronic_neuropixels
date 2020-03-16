% PLOT_PSYCHOMETRICS plot psychometric curves from the training and
% recording sessions
function [] = compare_psychometrics()
P = get_parameters;
T = readtable(P.performance_by_rat_path);
x = -30:30;

unique_rats = unique(T.rat);
for i = 1:numel(unique_rats)
    figure('Position', P.figure_position_psychometrics)
    set(gca, P.axes_properties_behavior{:}, ...
             'XLim', [min(x), max(x)], ...
             'Ylim', [0, 100]);
    plot(xlim, max(ylim)*[1,1], 'k--', 'LineWidth', 0.5)
    plot(xlim, 50*[1,1], 'k--', 'LineWidth', 0.5)
    plot(0*[1,1], ylim, 'k--', 'LineWidth', 0.5)
    idx1 = strcmp(T.rat, unique_rats{i}) & T.tethered;
    idx2 = strcmp(T.rat, unique_rats{i}) & ~T.tethered;
    y1 = eval_logistic4(x, T.gamma0(idx1), T.gamma1(idx1), T.sens(idx1), T.bias(idx1));
    y2 = eval_logistic4(x, T.gamma0(idx2), T.gamma1(idx2), T.sens(idx2), T.bias(idx2));
    plot(x,y1, 'k-', 'LineWidth', 1);
    plot(x,y2, 'k--', 'LineWidth', 1);
    xlabel('#R-#L clicks')
    ylabel('Percent chose right')
    title(unique_rats{i})
end