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
function [pval, T2] = anova_DV_shank(Cells, varargin)
parseobj = inputParser;
addParameter(parseobj, 'DV_bin_edges', [-10, -2, 0], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'EI_bin_edges', [1 577 768], ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
addParameter(parseobj, 'min_elec', 1, @(x) isscalar(x) && isnumeric(x));                                   
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P = get_parameters;

T = get_metrics_from_Cells(Cells, 'condition_on', {'DV', 'electrode_index'}, ...
                                  'DV_bin_edges', P_in.DV_bin_edges, ...
                                  'EI_bin_edges', P_in.EI_bin_edges, ...
                                  'x0', 0);
if strcmp(P_in.metric, 'Vpp')
    T.metric_norm = T.Vpp;
else
    T.metric_norm = T.(P_in.metric)./T.n_elec;
end
T = T(T.n_elec > P_in.min_elec, :);
%% Group
k = 0;
unique_id = unique(T.identifier);
for i = 1:numel(unique_id)
for i_DV = 1:numel(P_in.DV_bin_edges)-1
for i_EI = 1:numel(P_in.EI_bin_edges)-1
    idx = T.identifier==unique_id(i) & ...
          T.i_DV == i_DV & ...
          T.i_EI == i_EI;
    idx_early = idx & T.days_elapsed <= 7;
    idx_late = idx & T.days_elapsed > 14 & T.days_elapsed < 141;
    idx_late = idx & T.days_elapsed > 7;
    if sum(idx_early)==0 || sum(idx_late)==0
        continue
    end
    metric_early = max(T.metric_norm(idx_early)); 
    metric_late = mean(T.metric_norm(idx_late));
    metric_change = (metric_early - metric_late) / metric_early;
    k=k+1;
    T2.metric_change(k,1) = metric_change;
    T2.i_DV(k,1) = i_DV;
    T2.i_EI(k,1) = i_EI;
    T2.identifier(k,1) = unique_id(i);
end
end
end    
T2=struct2table(T2);
T2 = T2(~isnan(T2.metric_change),:);
%% Do the ANOVA
pval = anovan(T2.metric_change, {T2.i_EI, T2.i_DV}, 'continuous', [], 'display', 'off');