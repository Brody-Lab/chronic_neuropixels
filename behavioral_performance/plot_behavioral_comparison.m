% PLOT_BEHAVIORAL_COMPARISON compare a behavioral metric between training
% and tethered for recording
%
%=INPUT
%
%   data_type
%       A char array specifying the behavioral metric to be plotted
%
%=OPTIONAL INPUT, NAME-VALUE PAIRS
%
%   axes
%       The handle of the axes where the plot is made
%
%   label_n
%       Label the 
%
%   pval_pos
%       The relative [x,y] position from the axes's left, bottom corner. A
%       value of [0,0] is the lower left corner and a value of [1,1] is the
%       upper right corner. If it is empty, the pval is not shown.
%
%   samp_size_pos
%       The relative [x,y] position from the axes's left, bottom corner. A
%       value of [0,0] is the lower left corner and a value of [1,1] is the
%       upper right corner. If it is empty, the sample size is not shown.

function[]=plot_behavioral_comparison(data_type, varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'pval_pos', [0.7, 0.2], @(x) numel(x)==2||numel(x)==0);
addParameter(parseobj, 'samp_size_pos', [], @(x) numel(x)==2||numel(x)==0);
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P = get_parameters;
T = readtable(P.performance_by_rat_path);
if isempty(P_in.axes)
    figure('Position', P.figure_position_behavioral_comparison)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:}, ...
         'DataAspectRatio', [1,1,1])
if isfield(P.custom_axes_properties, data_type)
    set(gca, P.custom_axes_properties.(data_type){:})
end
unique_rats = unique(T.rat);
for i=1:numel(unique_rats)
    i_tethered = strcmp(T.rat, unique_rats{i}) & T.tethered;
    i_untether = strcmp(T.rat, unique_rats{i}) & ~T.tethered;
    data_tethered(i)=T.(data_type)(i_tethered);
    data_untether(i)=T.(data_type)(i_untether);
    if strcmp(data_type, 'bias')
        data_tethered(i)=abs(data_tethered(i));
        data_untether(i)=abs(data_untether(i));
    end
    plot(data_untether(i), data_tethered(i), 'o', 'Color', P.color_order(i,:));
end
diag_coords = [min([xlim,ylim]), max([xlim,ylim])];
plot(diag_coords,diag_coords, 'k','linewidth', 0.5);
xlabel('Untethered')
ylabel('Tethered for recording')
switch data_type
    case 'trials_done'
        title_text = 'Trials completed';
    case 'sens'
        title_text = 'Psychometric steepness';
    case 'abs_bias'
        title_text = '|Bias|';
    case 'lapse'
        title_text = 'Lapse rate';
    case 'prct_correct'
        title_text = 'Percent correct';
end
title(title_text)
% label pval
if ~isempty(P_in.pval_pos)
    pval = num2str(signrank(data_tethered, data_untether), '%.3f');
    ax_pos = get(gca, 'Position');
    x = P_in.pval_pos(1);
    y = P_in.pval_pos(2);
    annot_pos = [ax_pos(1)+ax_pos(3)*x, ...
                 ax_pos(2)+ax_pos(4)*y, ...
                 0, 0];
    annotation('textbox', annot_pos, 'String', ['p=' pval])
end
% label sample size
if ~isempty(P_in.samp_size_pos)
    n = numel(data_untether);
    ax_pos = get(gca, 'Position');
    x = P_in.samp_size_pos(1);
    y = P_in.samp_size_pos(2);
    annot_pos = [ax_pos(1)+ax_pos(3)*x, ...
                 ax_pos(2)+ax_pos(4)*y, ...
                 0, 0];
    annotation('textbox', annot_pos, 'String', ['N=' num2str(n)])
end