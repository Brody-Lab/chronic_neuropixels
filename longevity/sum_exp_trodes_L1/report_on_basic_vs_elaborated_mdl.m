function [prct_larger, pval] = report_on_basic_vs_elaborated_mdl(make_plot)
% REPORT_ON_BASIC_VS_ELABORATED_MDL
%
%=INPUT
%
%   make_plot
%       A logical indicating whether to make a plot
%
%=OUTPUT
%   
%   prct_larger
%       The percent increase in log-likelihood of the elaborated model
%       relative to the basic model
%
%   pval
%       probability that the difference in log-likelihood between the basic
%       and elaborated models is observed by chance

    add_folders_to_path
    P=get_parameters;
    load(P.sum_exp_trodes.orig_vs_elab_data_path, 'diff_boot', 'diff_obsv')
    
    prct_larger = (exp(diff_obsv)-1)*100;
    pval = sum(abs(diff_boot)>=abs(diff_obsv))/numel(diff_boot);
    pval = max(pval, 1/numel(diff_boot));
    
    if nargin < 1
        make_plot = false;
    end
    if make_plot
        figure
        set(gca, 'NextPlot', 'Add', ...
                 'fontsize', 14)
        hdl_boot = histogram(diff_boot);
        hdl_obsv = plot(diff_obsv*[1,1], ylim);
        xlabel('LL(elaborated) - LL(basic)')
        ylabel('# bootstrap samples')
        legend([hdl_boot, hdl_obsv], {'Boostrapped samples', 'Observed'}, 'location', 'best')
        title(['p \leq ' pval])
    end
end