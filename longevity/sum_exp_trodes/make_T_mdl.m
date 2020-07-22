% MAKE_T_MDL make a table of model variants with a unique subset of
% regressors
%
%=OPTIONAL INPUT
%
%    regressors
%       A string or cell array specifying the names of the regressors. The
%       regressors in the N1 term should start with "N1_" and those in the
%       change rate term should stat with "k_". A regressors whose name
%       "N1_xxx" is a replica of T_trode.(xxx). If xxx is "const," then a
%       vector of ones is created instead.
function T_mdl = make_T_mdl(varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'model_parameters', P.sum_exp_trodes.default_model_parameters, ...
        @(x) all(ismember(x, P.sum_exp_trodes.possible_model_parameters)))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    % number of regressors that are varied
    n_varying = sum(~contains(P_in.model_parameters, 'const'));
    inds = dec2bin(0:(2^n_varying-1));
    inds = inds=='1';

    k = 0;
    n_parameters = numel(P_in.model_parameters);
    for i = 1:n_parameters
        if contains(P_in.model_parameters{i}, 'const')
            T_mdl.(P_in.model_parameters{i}) = true(size(inds,1),1);
        else
            k = k + 1;
            T_mdl.(P_in.model_parameters{i}) = inds(:,k);
        end
    end
    T_mdl = struct2table(T_mdl);

    % remove model variants with both the categorical and the noncategorical
    % instance of the same model parameter
    idx_remove = false(size(T_mdl,1),1);
    for i = 1:n_parameters
        [i_regressor, ~, i_cat] = parse_regressor_name(P_in.model_parameters{i});
        for j = i+1:n_parameters
            [j_regressor, ~, j_cat] = parse_regressor_name(P_in.model_parameters{j});
            if strcmp(i_regressor, j_regressor) && (i_cat~=j_cat)
               idx_remove = idx_remove | all(T_mdl{:, [i,j]},2);
            end
        end
    end
    T_mdl = T_mdl(~idx_remove,:);
end
%% PARSE_REGRESSOR_NAME
%
% return the name, term in the equation where it appears, and whether it is
% categorical
%
%=INPUT
%   str
%       A char array of the model parameter name
%
%=OUTPUT
%
%   regressor
%       A char array specifying the regressor corresponding to
%       the model parameter
%
%   eq_term
%       A char array of the term in the equation where the regresssor
%       appears
%
%   is_cat
%       A logical indicating whether the regresssor is a categorical
%       variable
function [exp_factor, eq_term, is_cat] = parse_regressor_name(str)
    idx = regexp(str, '_*_');
    if numel(idx) > 1
        exp_factor =str(idx(1)+1:idx(2)-1);
        is_cat = true;
        eq_term = str(1:idx(1)-1);
    elseif numel(idx)==1
        exp_factor =str(idx(1)+1:end);
        is_cat = false;
        eq_term = str(1:idx(1)-1);
    elseif strcmp(str, 'alpha')
        exp_factor =str;
        is_cat = false;
        eq_term = str;
    else
        error('Error parsing regressor names')
    end
end