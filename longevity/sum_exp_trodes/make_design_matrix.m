% MAKE_DESIGN_MATRIX make a design matrix (table) from T_TRODE
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
%       A table that renames the experimental factors in T_trode according
%       to the optional input REGRESSORS
%
%=OPTIONAL INPUT
%
%    regressors
%       A string or cell array specifying the names of the regressors. The
%       regressors in the N1 term should start with "N1_" and those in the
%       change rate term should stat with "k_". A regressors whose name
%       "N1_xxx" is a replica of T_trode.(xxx). If xxx is "const," then a
%       vector of ones is created instead.
function T_dsgn = make_design_matrix(T_trode, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'regressors', P.sum_exp_trodes.regressors, @(x) all(ismember(x, P.longevity_metrics)))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
n = size(T_trode,1);
T_trode_var = T_trode.Properties.VariableNames;
for i = 1:numel(P_in.regressors)
    if contains(P_in.regressors{i}, 'const')
        T_dsgn.(P_in.regressors{i}) = ones(n,1);
    else
        idx_ = regexp(P.sum_exp_trodes.regressors{i}, '_');
        exp_fact_name = P.sum_exp_trodes.regressors{i}(idx_+1:end);
        idx = ismember(T_trode_var, exp_fact_name);
        if sum(idx)~=1
            error('Cannot find unique variable in T_TRODE corresponding to the regressor %s', ...
                   P.sum_exp_trodes.regressors{i})
        end
        T_dsgn.(P_in.regressors{i}) = T_trode{:,idx};
    end
end
T_dsgn = struct2table(T_dsgn);