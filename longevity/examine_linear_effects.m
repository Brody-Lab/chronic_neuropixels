% EXAMINE_LINEAR_EFFECTS test the effect of brain position and shank
% position on the number of units
%
%   Build a linear model with the regressors DAYS_ELAPSED, DV, and
%   AP, and SVP, and the interaction effects between DAYS_ELAPSED and
%   each other regressor.
function [mdl, T_trode] = examine_linear_effects(Cells, varargin)
parseobj = inputParser;
P = get_parameters;
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp', ...
                                                           'frac_single'}))
addParameter(parseobj, 'regressors', {'AP', 'DV', 'ML', 'SVP', 'SPA'})
addParameter(parseobj, 'x0', min(P.longevity_time_bin_centers), @(x) isscalar(x) && isnumeric(x))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
T_trode = struct;
for i =1:numel(Cells)
    if Cells{i}.days_since_surgery < P_in.x0
        continue
    end
    n_trode = sum(Cells{i}.electrodes.in_brain);
    T_trode.days_elapsed{i,1} = repmat(Cells{i}.days_since_surgery, n_trode, 1);
    T_trode.DV{i,1} = Cells{i}.electrodes.DV;
    T_trode.AP{i,1} = Cells{i}.electrodes.AP;
    T_trode.ML{i,1} = Cells{i}.electrodes.ML;
    T_trode.SVP{i,1} = ceil(Cells{i}.electrodes.index/2)*2/100; % um
    T_trode.SPA{i,1} = ones(n_trode,1)*(Cells{i}.shank_plane=="sagittal");
    for j = 1:n_trode
        is_nearby = abs(Cells{i}.electrode - Cells{i}.electrodes.index(j)) <=1; % two nearby sites
        T_trode.unit{i,1}(j,1) = sum(is_nearby);
        T_trode.single_unit{i,1}(j,1) = sum(is_nearby & ...
                                            Cells{i}.ks_good(:));
    end
end
T_trode = structfun(@cell2mat, T_trode, 'uni', 0);
T_trode = struct2table(T_trode);
%% Attempt a glmfit
T_norm = T_trode;
for i = 1:numel(P_in.regressors)
    r = T_norm.(P_in.regressors{i});
    T_norm.(P_in.regressors{i}) = (r-min(r))/(max(r)-min(r));
end
n = size(T_trode,1);
X=[];
% for i = 1:numel(P_in.regressors)
%     X = [X, T_norm.(P_in.regressors{i})];
% end
X = [X, -(T_norm.days_elapsed-1)];
% for i = 1:numel(P_in.regressors)
%     X = [X, T_norm.(P_in.regressors{i}).*-(T_norm.days_elapsed-1)];
% end
y = T_trode.(P_in.metric);
X = ones(size(y));
[b, dev, stats] = glmfit(X, y, 'poisson', 'link', 'log')
%% Compare observed and predicted
y_hat = glmval(b, X, 'log');
days_elapsed_unique = unique(T_trode.days_elapsed);
for i = 1:numel(days_elapsed_unique)
    idx = T_trode.days_elapsed==days_elapsed_unique(i);
    y_i =  T_trode.(P_in.metric)(idx);
    y_avg(i,1) = mean(y_i);
    y_sem(i,1) = sem(y_i);
    y_hat_i = y_hat(idx);
    y_hat_avg(i,1) = mean(y_hat_i);
    y_hat_sem(i,1) = sem(y_hat_i);
end
figure
hold on
errorbar(days_elapsed_unique, y_avg, y_sem)
errorbar(days_elapsed_unique, y_hat_avg, y_hat_sem)
set(gca, 'xscale', 'log')

figure
errorbar(y_avg, y_hat_avg, y_hat_sem,y_hat_sem, y_sem, y_sem, 'o')
set(gca, 'DataAspectRatio', [1,1,1])
hold on
refline(1)
%% Plot coefficients
figure
errorbar(1:numel(b), abs(b), stats.se*1.96)
set(gca, 'yscale', 'log')
set(gca, 'xlim', [0.5, numel(b)+0.5], ...
         'xtick', 1:numel(b), ...
         'xticklabel', [{'1'}, P_in.regressors, {'t'},  P_in.regressors])
%% OLD, trashed
keyboard
%% Fit a linear model
eqn = [P_in.metric ' ~'];
for i = 1:numel(P_in.regressors)
    if i == 1
        eqn = [eqn ' days_elapsed*' P_in.regressors{i}];
    else
        eqn = [eqn ' + days_elapsed*' P_in.regressors{i}];
    end
end
mdl = fitlm(T_trode, eqn);
%% Fit a linear model and then fit a glm
T_x0 = T_trode(T_trode.days_elapsed==P_in.x0,:);
eqn = [P_in.metric ' ~'];
X0 = [];
for i = 1:numel(P_in.regressors)
    if i == 1
        eqn = [eqn ' ' P_in.regressors{i}];
    else
        eqn = [eqn '+ ' P_in.regressors{i}];
    end
    X0 = [X0, T_norm.(P_in.regressors{i})];
end
mdl = fitlm(T_x0, eqn);
X0 = [ones(size(T_trode,1),1), X0];
A_hat = X0 * mdl.Coefficients.Estimate;
T_norm = T_trode;
T_norm.(P_in.metric) = T_norm.(P_in.metric)./A_hat;
X1 = -(T_norm.days_elapsed-P_in.x0) .* X0;
[b, dev, stats] = glmfit(X1, T_norm.(P_in.metric), 'poisson');

% b0 = [mdl.Coefficients.Estimate; b];
% X2 = [X0, -(T_norm.days_elapsed-P_in.x0)];
% modelfun = @(b, X)(X(:,1:4)*b(1:4) .* exp(X(:,5)./X(:,1:4)*b(5:8)));
% [beta,R,J,CovB,MSE,ErrorModelInfo]=nlinfit(X2, T_norm.(P_in.metric), modelfun, b0)
%%
figure
errorbar(1:numel(b), abs(b), stats.se*1.96)
set(gca, 'yscale', 'log')
y_hat = glmval(b, X1, 'log');
days_elapsed_unique = unique(T_norm.days_elapsed);
for i = 1:numel(days_elapsed_unique)
    idx = T_norm.days_elapsed==days_elapsed_unique(i);
    y_i =  T_norm.(P_in.metric)(idx);
    y_avg(i,1) = mean(y_i);
    y_sem(i,1) = sem(y_i);
    y_hat_i = y_hat(idx);
    y_hat_avg(i,1) = mean(y_hat_i);
    y_hat_sem(i,1) = sem(y_hat_i);
end
figure
hold on
errorbar(days_elapsed_unique, y_avg, y_sem)
errorbar(days_elapsed_unique, y_hat_avg, y_hat_sem)
set(gca, 'xscale', 'log')

figure
errorbar(y_avg, y_hat_avg, y_hat_sem,y_hat_sem, y_sem, y_sem, 'o')
set(gca, 'DataAspectRatio', [1,1,1])
hold on
refline(1)