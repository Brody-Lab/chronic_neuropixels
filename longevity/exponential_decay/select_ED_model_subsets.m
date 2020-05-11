% SELECT_ED_MODEL_SUBSETS select all model variants for the exponentia decay
% model
%
%=OUTPUT
%  
%   T_mdl
%       A table whose each row is a model variant and whose each column is
%       a regressor
function T_mdl = select_ED_model_subsets()
    P = get_parameters;
    n_regressors = numel(P.ED_trode_regressors);
    inds = dec2bin(0:(2^n_regressors-1));
    inds = inds== '1';
    
    % remove the models that contain both the continuous and categorical
    % variant of the same regressor
    idx_remove = false(size(inds,1),1);
    for i = 1:numel(P.ED_excluded_pairs)
        idx_remove = idx_remove | all(inds(:, P.ED_excluded_pairs{i}),2);
    end
    inds(idx_remove,:) = [];
    n_mdls = size(inds,1);
    
    % rearrange to have the 'y0' parameters be in front and add the
    % constant terms
    idx_y0 = contains(P.ED_trode_regressors, 'y0');
    inds = [true(n_mdls,1), inds(:,idx_y0), true(n_mdls,1), inds(:,~idx_y0)];
    var_names = ['y0', P.ED_trode_regressors(idx_y0), ...
                  'k', P.ED_trode_regressors(~idx_y0)];
    T_mdl = array2table(inds, 'variablenames', var_names);