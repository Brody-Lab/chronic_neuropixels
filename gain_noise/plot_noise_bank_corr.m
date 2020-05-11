% PLOT_NOISE_BANK_CORR plot the correlation across banks
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = plot_noise_bank_corr(varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj ,'example_sn',0);
addParameter(parseobj,'plot_mode','rsquare');
addParameter(parseobj,'max_z',Inf);
parse(parseobj, varargin{:});
P_in = parseobj.Results;
max_z = P_in.max_z;
analyze_gain_noise_data;

%% plot
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
    switch P_in.plot_mode
        case 'rsquare'
            data{i} = bank_rsquare_z{i};
        case 'linear'
            data{i}=bank_corr_linear_z{i};
        case 'rank'
            data{i} = bank_corr_rank_z{i};
    end
                
    if probe_sn{i} == P_in.example_sn
        hdl = plot(days_implanted{i}, data{i}, 'ro--', 'linewidth', 1,'linestyle','none');
    else
        hdl = plot(days_implanted{i}, data{i}, 'ko--', 'linewidth', 1,'linestyle','none');        
    end
    if days_implanted{i} == 0
        set(hdl, 'marker', '*');
        hdl_new = hdl;
    else
        hdl_explanted = hdl;
    end
end
xlabel('Cumulative days implanted')
ylabel('Across-bank R^2');
%legend([hdl_new, hdl_explanted], {'New', 'Explanted'}, 'Location', 'Best')
%% stats: Compare between fresh and explanted probes
fprintf('\nComparing the across-bank correlation on the electrodes between the new and explanated probes:')
days_implanted_latest = cellfun(@(x) x(end), days_implanted);
idx_new = days_implanted_latest == 0;
bank_0_data_new = cat(1,bank_0_noise_z{idx_new});
bank_1_data_new = cat(1,bank_1_noise_z{idx_new});
bank_0_data_exp = cat(1,bank_0_noise_z{~idx_new});
bank_1_data_exp = cat(1,bank_1_noise_z{~idx_new});
switch P_in.plot_mode
    case 'linear'
        func = @(x,y)corr(x,y);
        mode_string = 'Pearson r';
    case 'rank'
        func = @(x,y)corr(x,y,'Type','Spearman');
        mode_string = 'Spearman r';        
    case 'rsquare'
        func = @(x,y)rsquare(x,y);
        mode_string = 'R^2';        
end
boots_new = bootstrp(1000,func,bank_0_data_new,bank_1_data_new);
boots_exp = bootstrp(1000,func,bank_0_data_exp,bank_1_data_exp);
fprintf('\nUNIMPLANTED: %s=%0.3f, 95%% CI = [%0.3f %0.3f] (%i probes, n = %i electrodes) \nEXPLANTED: %s = %0.3f, 95%% CI = [%0.3f %0.3f] (%i explantations, n = %i electrodes)\n', ...
            mode_string,...
            func(bank_0_data_new,bank_1_data_new), ...
            prctile(boots_new,2.5),...
            prctile(boots_new,97.5),...
            sum(cellfun(@(x)~isnan(x),data(idx_new))), ...
            numel(bank_0_data_new), ...
            mode_string,...            
            func(bank_0_data_exp,bank_1_data_exp), ...
            prctile(boots_exp,2.5),...
            prctile(boots_exp,97.5),...          
            sum(cellfun(@(x)~isnan(x),data(~idx_new))), ...
            numel(bank_0_data_exp));