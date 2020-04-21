% PLOT_NOISE_BANK_CORR plot the correlation across banks
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = all_probes_bank_noise_scatter(varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj, 'new',false);
addParameter(parseobj,'explanted',false);
parse(parseobj, varargin{:});
P_in = parseobj.Results;

max_z = Inf;
analyze_gain_noise_data;
%% plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:});

if P_in.explanted
    h=scatter(cat(1,bank_0_noise_z{~idx_new}),cat(1,bank_1_noise_z{~idx_new}),'.');    
elseif P_in.new
    h=scatter(cat(1,bank_0_noise_z{idx_new}),cat(1,bank_1_noise_z{idx_new}),'.');        
else
    h=scatter(cat(1,bank_0_noise_z{:}),cat(1,bank_1_noise_z{:}),'.');        
end
set(gca,P.axes_properties{:},'xscale','linear','yscale','linear');
l=refline(1,0);
if P_in.explanted
    h=scatter(cat(1,bank_0_noise_z{~idx_new}),cat(1,bank_1_noise_z{~idx_new}),'o');    
elseif P_in.new
    h=scatter(cat(1,bank_0_noise_z{idx_new}),cat(1,bank_1_noise_z{idx_new}),'o');        
else
    h=scatter(cat(1,bank_0_noise_z{:}),cat(1,bank_1_noise_z{:}),'o');        
end
    
h.SizeData=20;
h.MarkerFaceColor='r';
h.MarkerEdgeColor='w';
h.LineWidth=0.1;
l.LineStyle =':';
xlabel('Bank 0 Noise (\muV_R_M_S, z-scored)');
ylabel('Bank 1 Noise (\muV_R_M_S, z-scored)');
%set(gca,'xtick',[-2 0 2 4],'ytick',[-2 0 2 4],'xticklabel',{'-2' '0' '2'  '\geq4'},'yticklabel',{ '-2' '0' '2'  '\geq4'},'xlim',[-2 4],'ylim',[-2 4]);

