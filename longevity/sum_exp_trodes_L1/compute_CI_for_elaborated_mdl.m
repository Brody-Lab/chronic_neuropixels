function [] = compute_CI_for_elaborated_mdl(varargin)
% COMPUTE_CI_FOR_ELABORATED_MDL compute the confidence intervals for the
% elaborated model
%
%=OPTIONAL INPUT
%
%   KFold
%       The number of folds for cross-validation
%
%   nboot
%       The number of bootstrap samples to draw. For each bootstrap sample,
%       each type of model is fitted. 

    P = get_parameters;
    load(P.Cells_path);
    parseobj = inputParser;
    addParameter(parseobj, 'KFold', 2, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
    addParameter(parseobj, 'nboot', 200, ...
        @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    T_trode = make_T_trode(Cells, 'model_parameters', P.sum_exp_trodes.possible_model_parameters, ...
                                  'normalize_regressors', true, ...
                                  'unit_distance', P.unit_distance, ...
                                  'x0', P.x0);
    n_trodes = size(T_trode,1);
    Selab = fit_SoE_L1(Cells, 'model_parameters', P.sum_exp_trodes.params_N1f_N1s, ...
                              'T_trode', T_trode);
    S.betas_obsv = Selab.T_betas;
    S.betas_boot = [];
    for i = 1:P_in.nboot
        fprintf('\n%i', i); tic

        idx = datasample(1:n_trodes, n_trodes);
        Selab = fit_SoE_L1(Cells, 'model_parameters', P.sum_exp_trodes.params_N1f_N1s, ...
                                  'T_trode', T_trode(idx,:));
        S.betas_boot = [S.betas_boot; Selab.T_betas];
        fprintf(' - took %0.f seconds', toc)
    end
    S.cil = quantile(S.betas_boot{:,:},0.025);
    S.ciu = quantile(S.betas_boot{:,:},0.975);
    S.lower = S.betas_obsv{:,:} - S.cil;
    S.upper = S.ciu - S.betas_obsv{:,:};
    
    save(P.sum_exp_trodes.elab_mdl_CI_path, 'S')
end