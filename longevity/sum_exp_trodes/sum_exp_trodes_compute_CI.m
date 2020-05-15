% SUM_EXP_TRODES_COMPUTE_CI compute the confidence interval for the
% coefficients estimated in the electrodes-sum-of-exponentials model
%
% Because the model selection is repeated multiple times, each time with a
% randomly selected half of the data, the subset associated with the median
% out-of-sample MSE is selected for calculating confidence intervals.
%
%=INPUT
%
%   S
%       A structure made using SELECT_EXP_MDL
%
%=OUTPUT
%
%   S
%       The same strcture but with the confidence intervals in S.T_res
%
%
%=OPTIONAL INPUT
%
%   n_boot
%       Number of bootstrap redraws
%
%   i_mdl
%       The index of the models whose parameters are to be estimated
%
%   iteration
%       Which iteration from which the samples are drawn. If this is empty,
%       then a random iteration is used.
function S = sum_exp_trodes_compute_CI(S, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'n_boot',P.exp_decay_n_boots, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
    addParameter(parseobj, 'i_mdl',1, ...
        @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results; 
    
    % add field
    n_regressors = size(S.T_res.b{1},1);
    n_mdls = size(S.T_res,1);
    for f = {'b_med', 'cil', 'ciu'}
        f = f{:};
        if ~isvar(S.T_res, f)
            S.T_res.(f) = nan(n_mdls, n_regressors);
        end
    end    
    if ~isvar(S.T_res, 'bootcoeff')
        S.T_res.bootcoef = cell(n_mdls,1);
    end
    if ~isvar(S.T_res, 'i_coef_est')
        S.T_res.i_coef_est = nan(n_mdls,1);
    end
    
    % for each model
    for m = 1:numel(P_in.i_mdl)
        ind_m = P_in.i_mdl(m);
        % use the subset of trials with median MSE
        [~,i] = min(S.T_res.MSE(ind_m,:)-median(S.T_res.MSE(ind_m,:)));
        bootcoef = nan(P_in.n_boot, n_regressors);
        i_trodes = ~S.trodes_used(:,i);
        n_trodes = sum(i_trodes);
        
        y =  S.T_trode.(S.P_in.metric)(i_trodes,:);
        t = S.T_trode.days_since_init(i_trodes);
        [XN1, Xk] = partition_T_dsgn(S.T_dsgn, S.T_mdl, ind_m, 'i_trodes', i_trodes);            
        parfor j = 1:P_in.n_boot
            bootidx = datasample(1:n_trodes, n_trodes);
            betas = fit_mdl_sum_exp_trodes(XN1(bootidx,:), ...
                             Xk(bootidx,:),t(bootidx,:), y(bootidx,:), 'noise', S.P_in.noise);

            b_tmp = nan(1, n_regressors);                
            b_tmp([true, true, true, S.T_mdl{ind_m,:}]) = betas;
            bootcoef(j,:) = b_tmp;
            fprintf('\n model %i boot %i',ind_m, j)
        end
        b = S.T_res.b{ind_m}(:, i);
        b=b(:)';
        S.T_res.b_med(ind_m,:) = b;
        S.T_res.cil(ind_m,:) =  quantile(bootcoef, 0.025);
        S.T_res.ciu(ind_m,:) =  quantile(bootcoef, 0.975);
        S.T_res.bootcoef{ind_m,1} = bootcoef;
        S.T_res.i_coef_est(ind_m,1) = i;
    end
end