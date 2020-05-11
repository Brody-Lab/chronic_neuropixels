% PARTITION_DESIGN_MATRIX partition a design matrix into separate matrices
% in each term of the exponential function
%
%=INPUT
%
%   T_dsgn
%       The design matrix (table) made by MAKE_T_DSGN
%
%   T_mdl
%       The table of model variants made by MAKE_T_MDL
function [XN1, Xk] = partition_T_dsgn(T_dsgn, T_mdl, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'i_mdl', true(size(T_mdl,1),1), ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {}))
addParameter(parseobj, 'i_trodes', true(size(T_dsgn,1),1), ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
T_dsgn = T_dsgn(i_trodes,:);
i_N1 = contains(T_mdl.Properties.VariableNames, 'N1_');
i_k = contains(T_mdl.Properties.VariableNames, 'k_');
XN1 = T_dsgn{P_in.i_trodes, T_mdl(i_N1)};
Xk  = T_dsgn{P_in.i_trodes, T_mdl(i_k)};