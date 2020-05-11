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
% The unit counts for each region are normalized by the average on the
% first day after the surgery. This normalization removes one fewer
% parameter to estimate. 

%% Parameters
re_bootstrap = false;
%% Get files
P = get_parameters;
if re_bootstrap
    if exist('Cells', 'var')
        assemble_exp_decay_data(Cells)
    else
        assemble_exp_decay_data()
    end
end
load(P.exp_decay_data_path)
%%
make_exp_decay_figure(T_ed, 'metric', 'unit')
make_exp_decay_figure(T_ed, 'metric', 'single_unit')