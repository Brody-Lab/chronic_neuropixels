% Make panels N and O of figure 2

if ~exist('S', 'var')
    P=get_parameters;
    load(P.sum_exp_trodes);
end
sum_exp_trodes_plot_coef(S)