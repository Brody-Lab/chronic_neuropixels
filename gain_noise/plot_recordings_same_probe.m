% PLOT_RECORDINGS_SAME_PROBE plot quality metrics from the same probe that
% was implanted in three animals
%
%=INPUT
%   
%   Cells
%       A cell array that was created by COLLECT_CELL_FILES
%
%=OPTIONAL INPUT
%
%   axes
%       The Axes object in which the plot is made. If empty (default), then
%       a new figure is created.
%
%   metric
%       A char array specifying the metric to be ploted
%
%   legend
%       A logical scalar specifying whether to show the legend
function [] = plot_recordings_same_probe(Cells, varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'metric', 'n_units', ...
    @(x) validateattributes(x, {'string', 'cell', 'char'}, {}))
addParameter(parseobj, 'legend', true, @(x) isscalar(x) && islogical(x));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
rat_name = ["T212"; "T224"; "T249"];
unique_bank = 0;
P = get_parameters;
T = make_table_from_Cells(Cells);
if isempty(P_in.axes)
    figure('Position', P.figure_position_longevity);
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:})
set(gca, P.custom_axes_properties.longevity{:});
for i = 1:numel(rat_name)
    idx=T.rat == rat_name(i) & T.unique_bank == unique_bank;
    % if there are multiple probes per animal, just use the one with the
    % largest number of recordings
    Ct = groupcounts(T(idx,:), 'probe_serial'); 
    idx = idx & T.probe_serial == Ct.probe_serial(end);    
%     T(idx,:) for debugging
    hdl(i) = plot(T.days_since_surgery(idx)+1, T.(P_in.metric)(idx), 'o-', ...
                    'Color', P.explant_color_order{i}, ...
                    'linewidth', 1);
    
    S.metric{i,1} = T.(P_in.metric)(idx);
    S.days_since_surgery{i,1}=T.days_since_surgery(idx);
    S.implant{i,1} = repmat(i, size(S.metric{i,1}));
end
set(gca, 'xlim', [2^-0.5, 2^7], ...
         'xtick', 2.^(0:2:6));
ylim([0,1].*ylim)
xlabel('Days since implant')
switch P_in.metric
    case 'n_units'
        ylabel('Units')
    case 'n_good_units'
        ylabel('Single units')
        set(gca,'ylim',[0 500]);
    case 'fr'
        ylabel('Event rate')
    case 'Vpp'
        ylabel({'Spike amplitude','(\muV_p_p)'});
end
if P_in.legend
    legend(hdl, {'Initial', 'Second', 'Third'}, 'location', 'northeast','FontSize',8,'color','w');
end
% Stats
S = structfun(@cell2mat, S, 'uni', 0);
S = struct2table(S);
pvals = anovan(S.metric, {S.implant}, 'display', 'off', 'continuous', [1]);
fprintf('\nANOVA on %s ', P_in.metric)
fprintf('\n   factor=the implant order: p = %0.3f', pvals(1))
% fprintf('\n   days since surgery: p = %0.3f', pvals(2))
