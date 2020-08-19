function X = make_design_matrix_for_each_term(T_dsgn, T_mdl, varargin)
% MAKE_DESIGN_MATRIX_FOR_EACH_TERM separate the design matrix for each term
%
%=INPUT
%
%   T_dsgn
%       The design matrix in the form of a table. Each column has a
%       variable name "x_y", where "x" indicates the term of the model and
%       "y" indicates the name of the regressor
%
%   T_mdl
%       A table specifying which model variables are used
%
%=OPTIONAL INPUT
%
%   i_mdl
%       The row in T_mdl that is selected
%
%   i_trodes
%       The rows in T_dsgn that are selected
%
%=OUTPUT
%
%   X
%       X is a structure whose fields are the names of the term 
parseobj = inputParser;
addParameter(parseobj, 'i_mdl', 1, ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {}))
addParameter(parseobj, 'i_trodes', true(size(T_dsgn,1),1), ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P_in.i_trodes = P_in.i_trodes(:);
X = struct;
parameter_names = T_mdl.Properties.VariableNames(T_mdl{P_in.i_mdl,:});
for i = 1:numel(parameter_names)
    parameter_type = get_parameter_type(parameter_names{i});
    regressor_name = get_regressor_name(parameter_names{i});
    if isfield(X, parameter_type)
        X.(parameter_type) = [X.(parameter_type), T_dsgn.(regressor_name)(P_in.i_trodes,:)];
    else
        X.(parameter_type) = T_dsgn.(regressor_name)(P_in.i_trodes,:);
    end
end
for param_type = {'N1', 'N1f', 'N1s', 'k'}
    if ~isfield(X, param_type{:})
        X.(param_type{:}) = [];
    end
end
assert(sum(T_mdl{P_in.i_mdl,:}) == sum(structfun(@(x) size(x,2), X)), 'Failed partition')