% MAKE_DESIGN_MATRIX make a design matrix (in a table form) from T_TRODE
%
%=INPUT
%
%   T_trode
%       Table made by MAKE_T_TRODE
%
%=OUTPUT
%
%   T_dsgn
%
%       A table that extracts the regressors from "T_trode" according to
%       "model_parameters". "T_dsgn" has the same number of rows as
%       "T_trode" and the same number of columns as the number of elements
%       of "model_parameters."
%
%=OPTIONAL INPUT
%
%   model_parameters
%       A char vector, string array, or cell array of char specifying the
%       model parameters. Each parameter must be a member of
%       P.possible_model_parameters. The parameters in the N1 term should
%       start with "N1_", those in the change rate term should start with
%       "k_", and those in the fraction of rapidly disppearing fraction
%       should start with "a_".
function T_dsgn = make_design_matrix(T_trode, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.default_model_parameters, ...
    @(x) all(ismember(x, P.sum_exp_trodes.possible_model_parameters)))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
for i = 1:numel(P_in.model_parameters)
    if contains(P_in.model_parameters{i}, 'const')
        T_dsgn.(P_in.model_parameters{i}) = ones(size(T_trode,1),1);
    else
        regressor_name = get_regressor_name(P_in.model_parameters{i});
        idx = ismember(T_trode.Properties.VariableNames, regressor_name);
        if sum(idx)~=1
            error('Cannot find unique variable in T_TRODE corresponding to the regressor %s', ...
                   P_in.model_parameters{i})
        end
        T_dsgn.(P_in.model_parameters{i}) = T_trode{:,idx};
    end
end
T_dsgn = struct2table(T_dsgn);