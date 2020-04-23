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
if nargin < 1
    load(P.Cells_path)
else
    Cells = varargin{1};
end
P = get_parameters;
T_ed = [];
for cond = {'', 'DV', 'AP', 'EI'};cond=cond{:};
    T = get_metrics_from_Cells(Cells, 'condition_on', cond);
    for metric = {'unit', 'single_unit'}; metric=metric{:};
        T_ed_i = make_ed_table(T, 'metric', metric);
        T_ed = [T_ed; T_ed_i];
    end
end
save(P.exp_decay_data_path, 'T_ed')