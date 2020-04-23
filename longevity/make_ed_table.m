% MAKE_ED_TABLE Make a data of exponental decay models
%
%=INPUT
%
%   T
%       A table of recording sessions separated by conditions made by
%       GET_METRICS_FROM_CELLS
%
%=OUTPUT
%
%   T_ed
%       A table of model fits, parameter estimates, confidence intervals
%       (95%), bootstrapped parameter estimates
%
%=OPTIONAL INPUT
%
%   fit_initial_value
%       A scalar logical specifiying whether the initial value of the
%       response variable (i.e., p(1)) will be fit or assumed to be the
%       mean of the response variables at x==x0.
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   n_boot
%       A scalar integer specifying the number of bootstrap draws
%
%   normalize_by_electrodes
%       A scalar logical specifying whether to normalize the metrics by the
%       number of electrodes
%
%   normalize_initial_values
%       A scalar logical specifying whether to normalize the metrics by the
%       value of the metric on the first day examined ("x0") 
%
%   x0
%       The initial days, by default 1 and specified in GET_PARAMETERS
%
function T_ed = make_ed_table(T, varargin)
P = get_parameters;
parseobj = inputParser;
addParameter(parseobj, 'fit_initial_value', false, ...
        @(x) isscalar(x)&&(x==0||x==1))
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
addParameter(parseobj, 'n_boot',P.exp_decay_n_boots, ...
        @(x) validateattributes(x, {'numeric'}, {'numel', 1}))
addParameter(parseobj, 'normalize_by_electrodes', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'normalize_initial_value', false, @(x) isscalar(x)&&(x==0||x==1))
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), ...
        @(x) validateattributes(x, {'numeric'}, {'scalar'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
if max(T.condition)>1
    P_in.normalize_by_electrodes = true;
end
if P_in.normalize_by_electrodes && ~strcmp(P_in.metric, 'Vpp')
    metric=T.(P_in.metric)./T.n_elec;
else
    metric=T.(P_in.metric);
end
condition_names = get_condition_names(T);
T_ed = struct;
for i = 1:max(T.condition)
    T_i = T(T.condition == i, :);    
    x = T_i.days_elapsed;
    y = T_i.(P_in.metric);
    if P_in.normalize_by_electrodes && ~strcmp(P_in.metric, 'Vpp')
        y = y./T_i.n_elec;
    end
    % name the response variable
    resp_var = P_in.metric;
    if P_in.normalize_initial_value
        resp_var = ['norm. ' resp_var];
        y = y/mean(y(x==P_in.x0));
    end
    if P_in.normalize_by_electrodes && ~strcmp(P_in.metric, 'Vpp')
        resp_var = [resp_var '/n_elec'];
    end
    T_ed.resp_var(i,1) = string(resp_var);
    % name the condition
    T_ed.cond_name(i,1) = condition_names(i);
    % fit 
    [p_hat, p_CI, p_boot] = fit_sum_2_exp_decay(x,y, ...
                                                'n_boot', P_in.n_boot, ...
                                                'x0', P_in.x0);
    tau_mdl_boot = compute_model_time_constant(p_boot, 'x0', P_in.x0);
    tau_mdl_boot = double(tau_mdl_boot);
    T_ed.tau_mdl(i,1) = median(tau_mdl_boot);
    % add data to the table of models
    T_ed.N0(i,1)    = p_hat(1);
    T_ed.alpha(i,1) = p_hat(2);
    T_ed.tau_f(i,1) = p_hat(3);
    T_ed.tau_s(i,1) = p_hat(4);
    T_ed.N0_CI(i,:)    = p_CI(:,1);
    T_ed.alpha_CI(i,:) = p_CI(:,2);
    T_ed.tau_f_CI(i,:) = p_CI(:,3);
    T_ed.tau_s_CI(i,:) = p_CI(:,4);
    T_ed.tau_mdl_CI(i,:) = quantile(tau_mdl_boot, [0.025, 0.975]);
    T_ed.N0_boot{i,1}    = p_boot(:,1);
    T_ed.alpha_boot{i,1} = p_boot(:,2);
    T_ed.tau_f_boot{i,1} = p_boot(:,3);
    T_ed.tau_s_boot{i,1} = p_boot(:,4);
    T_ed.tau_mdl_boot{i,1} = tau_mdl_boot;
    T_ed.metric(i,1) = string(P_in.metric);
    T_ed.normalize_by_electrodes(i,1) = P_in.normalize_by_electrodes;
    T_ed.fit_initial_value(i,1) = P_in.fit_initial_value;
    T_ed.T{i,1} = T;
    T_ed.condition(i,1) = i;
end
T_ed = struct2table(T_ed);