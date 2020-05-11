% FIT_EXP_TO_ELECTRODES fit an exponential decay function to the unit (or
% single_unit) count on each electrodes as a function of brain position or
% probe shank position of that electrode.
%
%=INPUT
%
%   Cells
%       The structure made by COLLECT_CELLS_FILE
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   n_boot
%       A scalar integer specifying the number of bootstrap draws. If
%       n_boot == 0, then the fits are based on the complete dataset
%
%   unit_distance
%       A scalar nonnegative, integer indicating the maximum distance (in
%       electrode index) between an electrode and the site where the peak
%       amplitude of a unit was observed for that unit to be counted for
%       that electrode. A value 0 indicates that a unit is counted for an
%       electrode only if that unit's peak response occurred on that
%       electrode.
%
%   regularization_alpha
%       The alpha for elastic net regularization in the interval (0,1]. A
%       value near 0 approaches ridge regularization, and 1 is lasso
%       regularization. If it's empty, not regularization will be
%       performed.
%
%   regressors
%       A char array specifying 
%
%   return_stats
%       A logical indicating whether to return the stats
%
%   x0
%       The first day after surgery to be examined.
function S = fit_exp_to_trodes(Cells, varargin)
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'hierarchical_bootstrap',false, @(x) isscalar(x) && (x==0 ||x==1))
addParameter(parseobj, 'unit_distance', 1, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
addParameter(parseobj, 'n_boot',P.exp_decay_n_boots, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
addParameter(parseobj, 'regularization_alpha',[],  ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>', 0, '<=', 1}))
addParameter(parseobj, 'regressors', P.exp_decay_regressors, ...
                @(x) all(ismember(x, P.exp_decay_regressors)))
addParameter(parseobj, 'whiten', [], @(x) all(ismember(x, {'PCA', 'ZCA', 'decorrelate'})))
addParameter(parseobj, 'return_stats', false, @(x) isscalar(x) && (x==0 ||x==1))
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), @(x) isscalar(x) && isnumeric(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
% Make T_trode
T_trode = make_T_trode(Cells, 'x0', P_in.x0, ...
                              'unit_distance', P_in.unit_distance);
% normalize the regressors
T_norm = T_trode; 
for i = 1:numel(P_in.regressors)
    r = T_norm.(P_in.regressors{i});
    T_norm.(P_in.regressors{i}) = (r-min(r))/(max(r)-min(r));
end
% whiten, as a control
if ~isempty(P_in.whiten)
    params = T_norm{:, P.exp_decay_regressors};
    S2 = cov(params);
    [V, D] = eig(S2);
    switch P_in.whiten
        case 'PCA'
            params = D^-0.5*V'*params';
        case 'ZCA'
            params = V*D^-0.5*V'*params';
        case 'decorrelate'
            params = V'*params';
    end
    T_norm{:, P.exp_decay_regressors} = params';
end
%% Get some variables
t = T_norm.days_elapsed-P_in.x0;
y = T_trode.(P_in.metric);
X = T_norm{:, P_in.regressors};
n = size(T_norm, 1);
ones_vec = ones(n,1);
%% FMINCON, additive effects on initial value
y0_hat = mean(y(T_norm.days_elapsed == P_in.x0));
n_reg = numel(P_in.regressors);
b0 = [y0_hat; zeros(n_reg,1); zeros(n_reg+1,1)];
X_ls = [ones_vec, X, ones_vec, X];
opts=optimset('fminunc');
opts.MaxFunEvals=1e6;
% opts.Display = 'off';
opts.TolFun = 1e-10;
b0=fminunc(@(b) ls_exp_additive(b, t, X_ls, y), b0, opts);
%% Bootstrap
k = 0;
S = struct;
betas = nan(P_in.n_boot, 2*(n_reg+1));
mse = nan(P_in.n_boot, 1);
exit_flags = nan(P_in.n_boot, 1);
parfor i = 1:P_in.n_boot
    boot_idx = datasample(1:n, n);
    [b,~,exit_flag] = fminunc(@(b) ls_exp_additive(b, t, X_ls(boot_idx,:), y(boot_idx)), b0, opts);
    fprintf('\n %i', i)
    betas(i,:)= b;
    mse(i,1) = ls_exp_additive(b,t,X_ls,y)/n;
    exit_flags(i,1) = exit_flag;
end
S.b = betas;
S.mse = mse;
S.exit_flag = exit_flags;
end

% y ~ exp(b{1}*x{1} + ... b{n}*x{n}) * ...
%     exp((b{n+1}*x{n+1} + ... b{2n}*x{2n})*t)
function Q = ls_exp_additive(b,t,X,y)
    y_hat = exp_additive(b,t,X);
    Q = sum((y - y_hat).^2)/size(X,1);
end

function y_hat = exp_additive(b,t,X)
    n_x = size(X,2)/2;
    y0 = X(:, 1:n_x)*b(1:n_x);
    k =  X(:, n_x+1:end)*b(n_x+1:end);
    y_hat = y0.*exp(k.*t);    
end