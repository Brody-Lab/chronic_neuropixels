% COMPARE_EXP_FIT_TO_TRODES Compare an exponential fit 
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
function S = compare_exp_fit_to_trodes(Cells, varargin)
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'KFold', 10, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
addParameter(parseobj, 'normalize_regressors',true, @(x) isscalar(x) && (x==0||x==1))
addParameter(parseobj, 'regressors', P.exp_decay_regressors, ...
                @(x) all(ismember(x, P.exp_decay_regressors)))
addParameter(parseobj, 'unit_distance', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), @(x) isscalar(x) && isnumeric(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
%% Make a table with information for each electrode recorded from each session
T_trode = make_T_trode(Cells, 'x0', P_in.x0, ...
                              'unit_distance', P_in.unit_distance);
for i = 1:numel(P_in.regressors)
    r = T_trode.(P_in.regressors{i});
    % need regressors to be positive for the model with multiplicative
    % effects of regressors on the initial value
    T_trode.(P_in.regressors{i}) = (r-min(r));
    if P_in.normalize_regressors
        T_trode.(P_in.regressors{i})=T_trode.(P_in.regressors{i})/(max(r)-min(r));
    end
end
%% name the basic some variables
t = T_trode.days_elapsed-P_in.x0;
y = T_trode.(P_in.metric);
X = T_trode{:, P_in.regressors};
n = size(T_trode, 1);
ones_vec = ones(n,1);
%% GLMFIT, multiplicative effect on initial value
% Make the design matrix and response variable
% Xt: the columns of the design matrix that depend on days_elapsed
X_glm =[X, [ones_vec, X] .* t];
%% FMINCON, additive effects on initial value
y0_hat = mean(y(T_trode.days_elapsed == P_in.x0));
n_reg = numel(P_in.regressors);
b0 = [y0_hat; zeros(n_reg,1); zeros(n_reg+1,1)];
X_ls = [ones_vec, X, ones_vec, X];

opts=optimset('fminunc');
opts.MaxFunEvals=1e6;
opts.TolFun = 1e-10;
%% Bootstrap
Boot_idx = bootstrp(P_in.n_boot, @(x)x, 1:size(T_trode,1)); 
k = 0;
S = struct;
for i = 1:P_in.n_boot
    fprintf('\n %i', i)
    boot_idx = Boot_idx(i,:)';
    
    b_glm = glmfit(X_glm(boot_idx,:), y(boot_idx), 'poisson');
    y_hat = glmval(b_glm,X_glm, 'log');
    k = k + 1;
    S.mse(k,1) = mean((y(boot_idx)-y_hat).^2);
    S.fit_type(k,1) = "glm_multi";
    S.boot(k,1) = i;
    
    b_ls_add = fminunc(@(b) ls_exp_additive(b, t, X_ls(boot_idx,:), y(boot_idx)), b0, opts);
    k = k + 1;
    S.mse(k,1) = ls_exp_additive(b_ls_add,t,X_ls,y)/n;
    S.fit_type(k,1) = "ls_addit";
    S.boot(k,1) = i;
    
    b_ls_mul = fminunc(@(b) ls_exp_multipl(b, t, X_ls(boot_idx,:), y(boot_idx)), b0, opts);
    k = k + 1;
    S.mse(k,1) = ls_exp_multipl(b_ls_mul, t, X_ls, y)/n;
    S.fit_type(k,1) = "ls_multi";
    S.boot(k,1) = i;
end
end
% y ~ exp(b{1}*x{1} + ... b{n}*x{n}) * ...
%     exp((b{n+1}*x{n+1} + ... b{2n}*x{2n})*t)
function Q = ls_exp_additive(b,t,X,y)
    y_hat = exp_additive(b,t,X);
    Q = sum((y - y_hat).^2);
end

function y_hat = exp_additive(b,t,X)
    n_x = size(X,2)/2;
    y0 = X(:, 1:n_x)*b(1:n_x);
    k =  X(:, n_x+1:end)*b(n_x+1:end);
    y_hat = y0.*exp(k.*t);    
end

% y ~ (b{1}*x{1} + ... b{n}*x{n}) * ...
%     exp((b{n+1}*x{n+1} + ... b{2n}*x{2n})*t)
function Q = ls_exp_multipl(b,t,X,y)
    y_hat = ex_multipl(b,t,X);
    Q = sum((y - y_hat).^2);
end

function y_hat = ex_multipl(b,t,X)
    n_x = size(X,2)/2;
    y0 = X(:, 1:n_x)*b(1:n_x);
    y0 = exp(y0);
    k =  X(:, n_x+1:end)*b(n_x+1:end);
    y_hat = y0.*exp(k.*t);    
end