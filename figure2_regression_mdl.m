% Make panels M,N,O,P of figure 2

add_folders_to_path

from_scratch = false;
if from_scratch
    sum_exp_trodes_run_mdl
else
    if ~exist('S', 'var')
        P = get_parameters;
        load(P.sum_exp_trodes.data_path);
    end
end

sum_exp_trodes_plot_coef(S, 'mdl_inds', 1:5)