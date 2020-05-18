% BOOTSTRAP_PVAL_TWO_SAMPLES P-value of the observed difference between two
% conditions in each parameter of the sum-of-exponentials model. 
%
% For each pair of groups, there are two samples x = {x_1, x_2, ..., x_m} and
% y = {y_1, y_2, ...y_n} of unit counts per session. Under the null hypothesis
% H_0: x = y, this calculates the probability that of the observed
% difference in each parameter theta of the sum-of-exponentials model
% between the two groups. A bootstrapping procedure is used.
%
% To estimate the null distribution, a new sample is defined z = {z_1, z_2,
% ..., z_{m+n}}. B bootstrap samples are drawn from this sample: x* =
% {x*_1, x*_2, ..., x*_m} and y* = {y*_1, y*_2, ..., y*_n}. For each
% bootstrap sample, the sum-of-exponentials model is fit to x* and y*. For
% each parameter theta, the test statistic t(x*,y*)_theta = theta_x* -
% theta_y*. The estimated probability value is #{t(x*,y*)_theta >=
% t(x,y)_theta. 
%
% If the table T provides information about more than two conditions, then
% the probability of each pair of conditions is computed.
%
%=INPUT
%
%   T
%       A table listing recording sessions created by
%       GET_METRICS_FROM_CELLS
%
%=OUTPUT
%
%   T_pval
%       A table showing the two-sample p-values
function T_pval = bootstrap_pval_two_samples(T, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'fit_initial_value', true, ...
        @(x) isscalar(x)&&(x==0||x==1))
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
addParameter(parseobj, 'n_boot',P.exp_decay_n_boots, ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 1}))
addParameter(parseobj, 'x0', P.x0, ...
        @(x) validateattributes(x, {'numeric'}, {'scalar'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
condition_names = get_condition_names(T);
k = 0;
regressor_names = {'N1', 'alpha', 'kf', 'ks'};
for i = 1:max(T.condition)-1
for j = i+1:max(T.condition)
    k =k + 1;
    T_pval.cond_i(k,1) = condition_names(i);
    T_pval.cond_j(k,1) = condition_names(j);
    t = T.days_elapsed;
    N = T.(P_in.metric);
    ti = t(T.condition==i);
    tj = t(T.condition==j);
    Ni = N(T.condition==i);
    Nj = N(T.condition==j);
    betas_obsv_i = fit_sum_2_exp_decay(ti,Ni, ...
                            'n_boot', 0, ...
                            'fit_initial_value', P_in.fit_initial_value, ...
                            'x0', P_in.x0);
    betas_obsv_j = fit_sum_2_exp_decay(tj,Nj, ...
                            'n_boot', 0, ...
                            'fit_initial_value', P_in.fit_initial_value, ...
                            'x0', P_in.x0);    
    % compute the bootstrap samples of the test-statistics
    m_init = sum(T.condition==i & T.days_elapsed == P.x0);
    n_init = sum(T.condition==j & T.days_elapsed == P.x0);
    m_subs = sum(T.condition==i & T.days_elapsed > P.x0);
    n_subs = sum(T.condition==j & T.days_elapsed > P.x0);
    
    t_init = t(T.days_elapsed == P.x0);
    t_subs = t(T.days_elapsed >  P.x0);
    N_init = N(T.days_elapsed == P.x0);
    N_subs = N(T.days_elapsed >  P.x0);
    
    betas_boot_i = nan(P_in.n_boot, 4);
    betas_boot_j = nan(P_in.n_boot, 4);
    parfor b = 1:P_in.n_boot
        rand_ind_init = datasample(1:m_init+n_init, m_init+n_init);
        rand_ind_subs = datasample(1:m_subs+n_subs, m_subs+n_subs);
        rand_ind_init_i = rand_ind_init(1:m_init);
        rand_ind_init_j = rand_ind_init(m_init+1:end);
        rand_ind_subs_i = rand_ind_subs(1:m_subs);
        rand_ind_subs_j = rand_ind_subs(m_subs+1:end);
        
        rand_ti = [t_init(rand_ind_init_i); t_subs(rand_ind_subs_i)];
        rand_tj = [t_init(rand_ind_init_j); t_subs(rand_ind_subs_j)];
        rand_Ni = [N_init(rand_ind_init_i); N_subs(rand_ind_subs_i)];
        rand_Nj = [N_init(rand_ind_init_j); N_subs(rand_ind_subs_j)];
        
        betas_boot_i(b,:) = fit_sum_2_exp_decay(rand_ti,rand_Ni, ...
                                    'n_boot', 0, ...
                                    'fit_initial_value', P_in.fit_initial_value, ...
                                    'x0', P_in.x0);
        betas_boot_j(b,:) = fit_sum_2_exp_decay(rand_tj,rand_Nj, ...
                                       'n_boot', 0, ...
                                       'fit_initial_value', P_in.fit_initial_value, ...
                                       'x0', P_in.x0);
    end
    % compare bootstrapped test-statistics to observed test-statistic
    for p = 1:numel(betas_obsv_i)
        if betas_obsv_i(p) > betas_obsv_j(p)
            ts_obsv = betas_obsv_i(p) - betas_obsv_j(p);
            tx_boot = betas_boot_i(:,p) - betas_boot_j(:,p);
        else
            ts_obsv = betas_obsv_j(p) - betas_obsv_i(p);
            tx_boot = betas_boot_j(:,p) - betas_boot_i(:,p);
        end            
        pval_1t = ones(1,P_in.n_boot)*(tx_boot >= ts_obsv)/P_in.n_boot;
        pval_2t = ones(1,P_in.n_boot)*(tx_boot >=  ts_obsv | ...
                                       tx_boot <= -ts_obsv)/P_in.n_boot;
        if pval_1t < eps
            pval_1t = 1/P_in.n_boot;
        end
        if pval_2t < eps
            pval_2t = 1/P_in.n_boot;
        end
        T_pval.(['pval_1t_' regressor_names{p}])(k,1) = pval_1t;
        T_pval.(['pval_2t_' regressor_names{p}])(k,1) = pval_2t;
    end
    
    % estimate the difference in the model time constant
    tau_mdl_obsv_i = compute_model_time_constant(betas_obsv_i, 'x0', P_in.x0);
    tau_mdl_obsv_j = compute_model_time_constant(betas_obsv_j, 'x0', P_in.x0);
    tau_mdl_boot_i = compute_model_time_constant(betas_boot_i, 'x0', P_in.x0);
    tau_mdl_boot_j = compute_model_time_constant(betas_boot_j, 'x0', P_in.x0);    
    if tau_mdl_obsv_i > tau_mdl_obsv_j
        ts_obsv = tau_mdl_obsv_i - tau_mdl_obsv_j;
        tx_boot = tau_mdl_boot_i - tau_mdl_boot_j;
    else
        ts_obsv = tau_mdl_obsv_j - tau_mdl_obsv_i;
        tx_boot = tau_mdl_boot_j - tau_mdl_boot_i;
    end    
    pval_1t = ones(1,P_in.n_boot)*(tx_boot >= ts_obsv)/P_in.n_boot;
    pval_2t = ones(1,P_in.n_boot)*(tx_boot >=  ts_obsv | ...
                                   tx_boot <= -ts_obsv)/P_in.n_boot;
    if pval_1t < eps
        pval_1t = 1/P_in.n_boot;
    end
    if pval_2t < eps
        pval_2t = 1/P_in.n_boot;
    end
    T_pval.('pval_1t_tau_mdl')(k,1) = pval_1t;
    T_pval.('pval_2t_tau_mdl')(k,1) = pval_2t;
    
    % save the estimated coefficients of the empirical samples and
    % bootstrap samples
    T_pval.betas_obsv_i(k,:) = betas_obsv_i;
    T_pval.betas_obsv_j(k,:) = betas_obsv_j;
    T_pval.betas_boot_i{k,1} = betas_boot_i;
    T_pval.betas_boot_j{k,1} = betas_boot_j;
    
    T_pval.tau_mdl_obsv_i(k,1) = tau_mdl_obsv_i;
    T_pval.tau_mdl_obsv_j(k,1) = tau_mdl_obsv_j;
    T_pval.tau_mdl_boot_i{k,1} = tau_mdl_boot_i;
    T_pval.tau_mdl_boot_j{k,1} = tau_mdl_boot_j;
    
    % save the input parameters
    for param_in = fieldnames(P_in)'; param_in=param_in{:};
        if ischar(P_in.(param_in))
            T_pval.(param_in)(k,1) = string(P_in.(param_in));
        else
            T_pval.(param_in)(k,1) = P_in.(param_in);
        end
    end
end
end
T_pval = struct2table(T_pval);