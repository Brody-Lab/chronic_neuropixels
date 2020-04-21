% PLOT_NOISE_BANK_CORR plot the correlation across banks
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = plot_across_bank_variance_explained(varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'new',false);
addParameter(parseobj,'explanted',false);
parse(parseobj, varargin{:});
P_in = parseobj.Results;

max_z=Inf;
analyze_gain_noise_data;

%% plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:}, ...
         'XLim', [0.5, 2.5], ...
         'Xtick', [1,2], ...
         'YLim', [0 1],...
         'XTicklabel', {'new', 'explanted'})

     
rsquare_boots_new = bootstrp(5000,@rsquare,cat(1,bank_1_noise_z{idx_new}),cat(1,bank_0_noise_z{idx_new}));
rsquare_boots_explanted = bootstrp(5000,@rsquare,cat(1,bank_1_noise_z{~idx_new}),cat(1,bank_0_noise_z{~idx_new}));

plot_data=[(rsquare_boots_new),(rsquare_boots_explanted)];

errorbar([1 2],mean(plot_data),mean(plot_data)-prctile(plot_data,2.5),-mean(plot_data)+prctile(plot_data,97.5));
   
ylabel('Across-bank variance explained');

