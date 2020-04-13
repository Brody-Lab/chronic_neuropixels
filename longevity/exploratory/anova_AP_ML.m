% anova_DV_shank dissociate the effect of DV
% position and shank position
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
%       of the electrode position. The second specifies the p-vale for the main effect of
%       the anatomical position
%
%   T2
%       A table summarizing the results on which the ANOVA was performed
%
%=OPTIONAL INPUT
%
%   DV_bin_edges
%       If CONDITION_ON includes 'DV', then this increasing vector
%       specifies the bin edges for binning DV positions
%
%   EI_bin_edges
%       If CONDITION_ON includes 'electode_index', then this increasing
%       vector specifies the bin edges for binning electrode indices. An
%       index of 1 indicates the electrode closest to the tip of the probe
%       shank, and an index of 960 indicates the electrode farthest from
%       the tip.
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   min_elec
%       Minimum number of electrodes for a recording/condition to be
%       used for the ANOVA
function [pval, T2] = anova_AP_ML(Cells, varargin)
parseobj = inputParser;
addParameter(parseobj, 'AP_bin_edges', [-10, -2, -1,  0], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
addParameter(parseobj, 'min_elec', 1, @(x) isscalar(x) && isnumeric(x));       
addParameter(parseobj, 'ML_bin_edges', [1 193, 385, 577, 768], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P = get_parameters;
T1 = get_metrics_from_Cells(Cells, 'condition_on', 'electrode_index', ...
                                   'EI_bin_edges', [1 193, 385, 577, 768]);
if strcmp(P_in.metric, 'Vpp')
    T1.metric_norm = T1.Vpp;
else
    T1.metric_norm = T1.(P_in.metric)./T1.n_elec;
end
T1 = T1(T1.n_elec > P_in.min_elec, :);
%% Group
k = 0;
unique_id = unique(T1.identifier);
unique_cond = unique(T1.condition); % this should be (1:n_cond)'
for i = 1:numel(unique_id)
for j = 1:numel(unique_cond)
    idx = T1.identifier==unique_id(i) & ...
          T1.condition == unique_cond(j);
    idx_early = idx & T1.days_elapsed <= 7;
    idx_late = idx & T1.days_elapsed > 28 & T1.days_elapsed < 141;
    if sum(idx_early)==0 || sum(idx_late)==0
        continue
    end
    metric_early = max(T1.metric_norm(idx_early)); 
    metric_late = mean(T1.metric_norm(idx_late));
    metric_change = (metric_early - metric_late) / metric_early;
    k=k+1;
    T2.metric_change(k,1) = metric_change;
    T2.AP(k,1) = unique(T1.trode_AP_avg(idx));
    T2.ML(k,1) = unique(T1.trode_ML_avg(idx));
    T2.identifier(k,1) = unique_id(i);
    T2.condition(k,1) = unique_cond(j);
end  
end
T2=struct2table(T2);
T2 = T2(~isnan(T2.metric_change),:);
%% Do the ANOVA
pval = anovan(T2.metric_change, {T2.AP, T2.ML}, 'continuous', [1,2], ...
                                                'display', 'off');