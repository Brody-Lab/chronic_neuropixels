% QUANTIFY_STABILITY_DEPENDENCE dissociate the effect of anatomical
% position and bank
%
%=INPUT
%
%   Cells
%       Structure made by COLLECT_CELL_FILES
%
%=OUTPUT
%
%   pval
%       The first of two element specifies the p-value for the main effect
%       of the bank. The second specifies the p-vale for the main effect of
%       the anatomical position
%
%   T2
%       A table summarizing the results on which the ANOVA was performed

function [pval, T2] = quantify_stability_dependence(Cells, varargin)
parseobj = inputParser;
addParameter(parseobj, 'anatom_bin_edges', [-10, -2, -1,  0], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'condition_on', 'DV', @(x) ismember(x, {'ML', ...
                                                             'AP', ...
                                                             'DV'}))
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P = get_parameters;
%% Get the appropriate metric
switch P_in.metric
    case 'single_unit'
        Metric = cellfun(@(x) x.ks_good,Cells,'uni',0);
        func = @sum;
    case 'unit'
        Metric = cellfun(@(x) ones(x.n_units,1), Cells, 'uni', 0);
        func = @sum;
    case 'event_rate'
        Metric = cellfun(@(x)x.fr,Cells,'uni', 0);
        func = @sum;
    case 'Vpp'
        Metric = cellfun(@(x) x.unitVppRaw,Cells,'uni', 0);
        func = @mean;
        P_in.normalize_by_electrodes = false;
end
y_label = P.text.(P_in.metric);
Metric = cellfun(@(x) x(:),Metric,'uni', 0);
%% Separate the metric from each recording
k = 0;
for i = 1:numel(Metric)
    for m = [0,1]
    for n = 1:numel(P_in.anatom_bin_edges)-1
        k=k+1;
        idx_cells = Cells{i}.bank == m & ...
                    Cells{i}.(P_in.condition_on) >= P_in.anatom_bin_edges(n) & ...
                  Cells{i}.(P_in.condition_on) <  P_in.anatom_bin_edges(n+1);
        idx_elec = Cells{i}.electrodes.bank == m & ...
                   Cells{i}.electrodes.(P_in.condition_on) >= P_in.anatom_bin_edges(n) & ...
                   Cells{i}.electrodes.(P_in.condition_on) <  P_in.anatom_bin_edges(n+1);
        T.identifier(k,1) = string(Cells{i}.identifier);
        T.days_elapsed(k,1) = Cells{i}.days_since_surgery;
        T.metric(k,1) = func(Metric{i}(idx_cells));
        T.n_elec(k,1) = sum(idx_elec);
        T.bank(k,1) = m;
        T.anatom_bin(k,1) = n;
    end
    end
end
if strcmp(P_in.metric, 'Vpp')
    T.metric_norm = T.metric;
else
    T.metric_norm = T.metric./T.n_elec;
end
T = struct2table(T);
T = T(T.n_elec > 0, :);
%% Group
k = 0;
unique_id = unique(T.identifier);
for i = 1:numel(unique_id)
for m = [0,1]
for n = 1:numel(P_in.anatom_bin_edges)-1
    idx = T.identifier==unique_id(i) & ...
          T.bank == m & ...
          T.anatom_bin == n;
    idx_early = idx & T.days_elapsed <= 7;
    idx_late = idx & T.days_elapsed > 28 & T.days_elapsed < 141;
    if sum(idx_early)==0 || sum(idx_late)==0
        continue
    end
    metric_early = max(T.metric_norm(idx_early)); 
    metric_late = mean(T.metric_norm(idx_late));
    metric_change = (metric_early - metric_late) / metric_early;
    k=k+1;
    T2.metric_change(k,1) = metric_change;
    T2.bank(k,1) = m;
    T2.anatom_bin(k,1) = n;
    T2.identifier(k,1) = unique_id(i);
end
end
end    
T2=struct2table(T2);
T2 = T2(~isnan(T2.metric_change),:);
%% Do the ANOVA
pval = anovan(T2.metric_change, {T2.bank, T2.anatom_bin}, 'continuous', [1,2], 'display', 'off');