function varargout = report_on_basic_vs_elaborated_mdl(make_plot)
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
    load(P.sum_exp_trodes.orig_vs_elab_data_path)
    
    prct_larger = (exp(diff_obsv)-1)*100;
    pval = sum(abs(diff_boot)>=abs(diff_obsv))/numel(diff_boot);
    pval = max(pval, 1/numel(diff_boot));
    
    fprintf('\nThe loglikelihood of the elaborated model is higher by %0.1f%%', ...
            prct_larger)
    fprintf(' (p = %0.3f)\n', pval)
    if nargout > 0
        varargout{1} = prct_larger;
    end
    if nargout > 1
        varargout{2} = pval;
    end
    
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
    %% BIC
    BIC_diff = median(Scomp.BIC_elab - Scomp.BIC_orig);
    if BIC_diff < 0
        fprintf('\n The elaborated SoER model has a lower BIC score (difference = %0.0f)', ...
                 BIC_diff)
    else
        fprintf('\n The basic SoER model has a lower BIC score (difference = %0.0f)', ...
                 -BIC_diff)
    end
end