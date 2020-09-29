function Scomp = compare_orig_elab_models(varargin)
% COMPARE_ORIG_ELAB_MODELS compare the out-of-sample log-likelihood of the
% original and the elaborated model. The original model was sparsified
% using L0 regularization. The elaborated model was sparsified using an L1
% regularization.
%
%=OPTIONAL INPUT
%
%   KFold
%       The number of folds for cross-validation
%
%   nboot
%       The number of bootstrap samples to draw. For each bootstrap sample,
%       each type of model is fitted. 
%
%=OUTPUT
%
%   Scomp
%       A structure containing information on the comparison "LL_elab" and
%       "LL_orig"
    P = get_parameters;
    load(P.Cells_path);
    parseobj = inputParser;
    addParameter(parseobj, 'KFold', 2, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'nboot', 20, ...
        @(x) all(ismember(x, P.sum_exp_trodes.possible_model_parameters)))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    T_trode = make_T_trode(Cells, 'model_parameters', P.sum_exp_trodes.possible_model_parameters, ...
                                  'normalize_regressors', true, ...
                                  'unit_distance', P.unit_distance, ...
                                  'x0', P.x0);
    n_trodes = size(T_trode,1);
    
    Scomp.i_opt_lambda_elab = nan(P_in.nboot,1);
    Scomp.LL_elab = nan(P_in.nboot,1);
    Scomp.i_opt_lambda_orig = nan(P_in.nboot,1);
    Scomp.LL_orig = nan(P_in.nboot,1);
    
    for i = 1:P_in.nboot
        fprintf('\n%i', i); tic
            
        idx = datasample(1:n_trodes, n_trodes);
        i_est_lambda = ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
        cvp = cvpartition(sum(i_est_lambda), 'KFold', 2);

        Selab = fit_SoE_L1(Cells, 'model_parameters', P.sum_exp_trodes.params_N1f_N1s, ...
                                  'T_trode', T_trode(idx,:), ...
                                  'i_est_lambda', i_est_lambda, ...
                                  'cvp', cvp, ...
                                  'infer_betas', false);
        Sorig = fit_SoE_L1(Cells, 'model_parameters', P.sum_exp_trodes.params_orig, ...
                                  'num_lambda', 1, ...
                                  'T_trode', T_trode(idx,:), ...
                                  'i_est_lambda', i_est_lambda, ...
                                  'cvp', cvp, ...
                                  'infer_betas', false);
        Scomp.i_opt_lambda_elab(i,1) = Selab.i_opt_lambda;
        Scomp.LL_elab(i,1) = Selab.T_lambda.LL_per_trode(Selab.i_opt_lambda);
        Scomp.i_opt_lambda_orig(i,1) = Sorig.i_opt_lambda;
        Scomp.LL_orig(i,1) = Sorig.T_lambda.LL_per_trode(Sorig.i_opt_lambda);
        Scomp.BIC_elab(i,1) = Selab.BIC;
        Scomp.BIC_orig(i,1) = Sorig.BIC;
        fprintf(' - took %0.f seconds', toc)
    end
    %% Do the second set of bootstrap
    null_LL_elab = Scomp.LL_elab - mean(Scomp.LL_elab);
    null_LL_orig = Scomp.LL_orig - mean(Scomp.LL_orig);
    diff_boot = bootstrp(1000, @(x,y) mean(x)-mean(y), null_LL_elab, null_LL_orig);
    
    % calculated the observed difference
    i_est_lambda = ismember(1:n_trodes, randperm(n_trodes, round(n_trodes/2)));
    cvp = cvpartition(sum(i_est_lambda), 'KFold', 2);
    Selab = fit_SoE_L1(Cells, 'model_parameters', P.sum_exp_trodes.params_N1f_N1s, ...
                                  'T_trode', T_trode, ...
                                  'i_est_lambda', i_est_lambda, ...
                                  'cvp', cvp);
    Sorig = fit_SoE_L1(Cells, 'model_parameters', P.sum_exp_trodes.params_orig, ...
                              'num_lambda', 1, ...
                              'T_trode', T_trode, ...
                              'i_est_lambda', i_est_lambda, ...
                              'cvp', cvp);
    diff_obsv = Selab.T_lambda.LL_per_trode(Selab.i_opt_lambda) - ...
                Sorig.T_lambda.LL_per_trode(Sorig.i_opt_lambda);
    %% Save results
    save(P.sum_exp_trodes.orig_vs_elab_data_path, 'Scomp', 'diff_boot', 'diff_obsv')
end