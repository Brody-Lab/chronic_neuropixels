% PLOT_GAIN_NOISE_BROKEN_FRAC The fraction of electrodes that are above a
% noise threshold
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = plot_gain_noise_broken_frac(varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
analyze_gain_noise_data;
%% fraction of electrodes with RMS noise above threshold
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:}, ...
         'XLim', [-10, 200], ...
         'Xtick', [0, 100, 200], ...
         'XTicklabel', {'0', '100', '\geq200'})
for i =1:n_probes
    if unique_probes(i)==17131312432
        hdl = plot(days_implanted{i}, frac_noisy{i}*100,  'linewidth', 1,'markersize',6,'color','k');        
        for z=1:length(days_implanted{i})
            scatter(days_implanted{i}(z),frac_noisy{i}(z)*100,40,'markeredgecolor',P.explant_color_order{z},'markerfacecolor','none','linewidth',1);
        end
    else
        hdl = plot(days_implanted{i}, frac_noisy{i}*100, 'ko--', 'linewidth', 1,'markersize',6);
    end
    if days_implanted{i} == 0
        set(hdl, 'marker', '*');
    end
end
set(gca, 'Ylim', ylim.*[0,1])
xlabel('Cumulative days implanted')
ylabel({'Percent >',[num2str(P.noise_threshold_uV) ' \muV_R_M_S']});
%% stats
idx_new = cell2mat(days_implanted) == 0;
noise_uV_all=vertcat(noise_uV{:});
noise_uV_new = cell2mat(noise_uV_all(idx_new));
noise_uV_exp = cell2mat(noise_uV_all(~idx_new));
x.new(1,1) = sum(noise_uV_new<P.noise_threshold_uV);
x.new(2,1) = sum(noise_uV_new>=P.noise_threshold_uV);
x.exp(1,1) = sum(noise_uV_exp<P.noise_threshold_uV);
x.exp(2,1) = sum(noise_uV_exp>=P.noise_threshold_uV);
x=struct2table(x);
[~,pval]=fishertest(x);
[~,ci_new] = binofit(x.new(2), sum(x.new));
[~,ci_exp] = binofit(x.exp(2), sum(x.exp));
fprintf('\nComparing the fraction of electrodes with a RMS noise > 20 uV between new and explanted probes:')
fprintf('\n   new: %0.2f,95%%CI: [%0.2f, %0.2f]', x.new(2)/sum(x.new)*100, ci_new(1)*100, ci_new(2)*100)
fprintf('\n   exp: %0.2f 95%%CI: [%0.2f, %0.2f]', x.exp(2)/sum(x.exp)*100, ci_exp(1)*100, ci_exp(2)*100)
fprintf('\n   (using the Fisher'' exact test):')
fprintf('\n   p = %f', pval)
