% COMPARE_MODEL_TYPES_FOR_TRODES Compare whether a linear model or an
% exponential model fit the unit count for each electrode 
%
%=INPUT
%
%   Cells
%       The structure made by COLLECT_CELLS_FILE and POSTPROCESS_CELLS
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   KFold
%       Number of cross-validation folds
%
%   MCReps
%       Number of Monte Carlo repetitions
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   unit_distance
%       A scalar nonnegative, integer indicating the maximum distance (in
%       electrode index) between an electrode and the site where the peak
%       amplitude of a unit was observed for that unit to be counted for
%       that electrode. A value 0 indicates that a unit is counted for an
%       electrode only if that unit's peak response occurred on that
%       electrode.
%
%   x0
%       The first day after surgery to be examined.
%
%=OUTPUT
%
%   S
%       A table with the fields "mse" and "fit_type"
function S = compare_model_types_for_trodes(Cells, varargin)
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'KFold', 10, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
addParameter(parseobj, 'MCReps', 10, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonzero'}))
addParameter(parseobj, 'metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
addParameter(parseobj, 'unit_distance', 0, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'integer', 'nonnegative'}))
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), @(x) isscalar(x) && isnumeric(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
T_trode = make_T_trode(Cells, 'x0', P_in.x0, ...
                              'unit_distance', P_in.unit_distance);
t = T_trode.days_elapsed-P_in.x0;
y = T_trode.(P_in.metric);
n = size(T_trode, 1);
X = [ones(n,1), t];
S.mse = [];
S.fit_type = [];
for i = 1:P_in.MCReps
    mse = crossval('mse', X, y, 'Predfun', @fp_lm, ...
                    'KFold', P_in.KFold);
    S.mse = [S.mse; mse];
    S.fit_type = [S.fit_type; "lm"];

    mse = crossval('mse', X, y, 'Predfun', @fp_poisson_glm, 'KFold', ...
                        P_in.KFold);
    S.mse = [S.mse; mse];
    S.fit_type = [S.fit_type; "poisson_glm"];

    mse = crossval('mse', X, y, 'Predfun', @fp_poisson_glm_fminunc,...
                    'KFold', P_in.KFold);
    S.mse = [S.mse; mse];
    S.fit_type = [S.fit_type; "poisson_glm_fminunc"];
end
S = struct2table(S);
end
%% fit a linear regression and predict values
function yfit = fp_lm(Xtrain, ytrain, Xtest)
    b = Xtrain\ytrain;
    yfit = Xtest*b;
end
%% fit a exponential-linear glm using glmfit and predict values
function yfit = fp_poisson_glm(Xtrain, ytrain, Xtest)
    b = glmfit(Xtrain, ytrain, 'poisson', 'constant', 'off');
    yfit = glmval(b, Xtest, 'log', 'constant', 'off');
end
%% fit an exponential function using fminunc and predict values
function yfit = fp_poisson_glm_fminunc(Xtrain, ytrain, Xtest)
    opts=optimset('fminunc');
    opts.MaxFunEvals=1e6;
    opts.TolFun = 1e-10;
    b = fminunc(@(b) ls_exp(b, Xtrain, ytrain), [1,-1], opts);
    yfit = eval_exp(b, Xtest);
end
%% Caclulate the sum of least squares
function Q = ls_exp(b,X,y)
    y_hat = eval_exp(b,X);
    Q = sum((y - y_hat).^2);
end
%% Evaluate a exponential function
function yfit = eval_exp(b,X)
    yfit = b(1)*X(:,1).*exp(b(2)*X(:,2));   
end