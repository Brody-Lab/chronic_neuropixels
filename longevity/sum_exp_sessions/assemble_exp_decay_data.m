% ASSEMBLE_EXP_DECAY_DATA assemble data for models of exponential decay
%
% The data are saved as a table at the location specified in
% GET_PARAMETERS, by the field "exp_decay_data_path"
%
% A table is created by MAKE_ED_TABLE for each condition, and the tables
% from all conditions are stacked. 
%
%=INPUT
%   
%   Cells
%       A structure made using COLLECT_CELLs_FILES and POSTPROCESS_CELLS
%
function [] = assemble_exp_decay_data(varargin)
P = get_parameters;
if nargin < 1
    load(P.Cells_path)
else
    Cells = varargin{1};
end
T_ed = [];
T_pval = []; 
for cond = {'', 'DV', 'AP'};cond=cond{:};
    T = get_metrics_from_Cells(Cells, 'condition_on', cond);
    for metric = {'unit', 'single_unit'}; metric=metric{:};
        T_ed_i = make_ed_table(T, 'metric', metric, ...
                                  'fit_initial_value', true);
        T_ed = [T_ed; T_ed_i];
        % If there are two condition associated with the samples X and Y,
        % a p-val is computed for each parameter theta that theta_X =
        % theta_Y under the null hypothesis H_0: X = Y. The p-val is
        % computed using bootstrapping.
        if max(T.condition)==2
            T_pval_i = bootstrap_pval_two_samples(T, 'metric', metric, ...
                                                     'fit_initial_value', true);
            T_pval=[T_pval;T_pval_i];
        end
    end
end
save(P.sum_exp_data_path, 'T_ed', 'T_pval')