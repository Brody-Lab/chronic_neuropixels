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
function [S, X, y] = fit_exp_to_electrodes(Cells, varargin)
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
T_trode = make_T_trode(Cells, 'x0', P_in.x0);
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
% Make the design matrix and response variable
% X0: the columns of the design matrix that do not depend on days_elapsed
X0 = T_norm{:, P_in.regressors};
% Xt: the columns of the design matrix that depend on days_elapsed
Xt =[ones(size(T_norm,1),1), T_norm{:, P_in.regressors}] .* (T_norm.days_elapsed-P_in.x0);
X = [X0, Xt];
y = T_trode.(P_in.metric);
% Bootstrap
if P_in.hierarchical_bootstrap
    unique_Cells_index = unique(T_norm.Cells_index);
    i = 0;
    while i<P_in.n_boot
        boot_idx = [];
        boot_sess = datasample(unique_Cells_index, numel(unique_Cells_index));
        for j = 1:numel(boot_sess)
            idx = T_norm.Cells_index == boot_sess(j);
            inds = datasample(find(idx), sum(idx));
            boot_idx = [boot_idx; inds];
        end
        lastwarn('');
        FitInfo = struct;
        if isempty(P_in.regularization_alpha)
            [b, dev, stats] = glmfit(X(boot_idx,:), y(boot_idx), 'poisson');
            FitInfo.dev=dev;
            FitInfo.stats=stats;
        else
            [b, FitInfo] = lassoglm(X(boot_idx,:), y(boot_idx), 'poisson', ...
                                    'Alpha', 0.01);
        end
        if isempty(lastwarn)
            i = i + 1;
        else
            continue
        end
        S(i).b = b;
        if P_in.return_stats
            S(i).FitInfo = FitInfo;
        end
        fprintf('\n%i', i)
    end
else
    if P_in.n_boot < 1 % or not
        Boot_idx = 1:size(T_norm,1);
    else
        Boot_idx = bootstrp(P_in.n_boot, @(x)x, 1:size(T_norm,1)); 
    end
    for i = 1:max(P_in.n_boot,1)
        boot_idx = Boot_idx(i,:)';
        FitInfo = struct;
        if isempty(P_in.regularization_alpha)
            [b, dev, stats] = glmfit(X(boot_idx,:), y(boot_idx), 'poisson');
            FitInfo.dev=dev;
            FitInfo.stats=stats;
        else
            [b, FitInfo] = lassoglm(X(boot_idx,:), y(boot_idx), 'poisson', ...
                                    'Alpha', P_in.regularization_alpha);
        end
        S(i).b = b;
        if P_in.return_stats
            S(i).FitInfo = FitInfo;
        end
    end
end