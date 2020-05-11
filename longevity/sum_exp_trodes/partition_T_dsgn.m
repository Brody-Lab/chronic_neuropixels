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
%
%   i_mdl
%       The model index in T_mdl
%
%=OUTPUT
%
%   XN1
%       The design matrix for the N1 term
%
%   Xk
%       The design matrix for the change rate (kappa) term
%
%=OPTIONAL INPUT
%
%   i_trodes
%       The electrodes that are being used
%
%=EXAMPLE
%

function [XN1, Xk] = partition_T_dsgn(T_dsgn, T_mdl, i_mdl, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'i_trodes', true(size(T_dsgn,1),1), ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P_in.i_trodes = P_in.i_trodes(:);
i_N1 = contains(T_mdl.Properties.VariableNames, 'N1_') & T_mdl{i_mdl, :};
i_k = contains(T_mdl.Properties.VariableNames, 'k_') & T_mdl{i_mdl, :};
XN1 = T_dsgn{P_in.i_trodes, i_N1};
Xk  = T_dsgn{P_in.i_trodes, i_k};
assert(sum(T_mdl{i_mdl,:}) == size(XN1,2) + size(Xk,2), 'Failed partition')