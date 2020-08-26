% FIGURE2_SUPP_SUM_EXP make the plots for a Figure 4--supplement that
% provides details about the 
%
% The model:
%
%   n_units ~ alpha*exp(-(x-x0)/tau1) + (1-alpha)*exp(-(x-x0)/tau2) + epsilon_i
%
% The parameter "alpha" is the fraction of units that shows more rapid loss
% and is between 0 and 1. The parameter "tau1" is the time constant
% associated with the more rapid degradation, and "tau2" is associated with
% the slower time constant. The parameter "x0" is the earliest day since
% the surgery when we include data for analysis, and is 1. Finally,
% "epsilon_i" is i.i.d. Gaussians with the same variance. I have not
% determined whether a heteroskedastic model might change the conclusions.
%
% The unit counts for each region are normalized by the average on the
% first day after the surgery. This normalization removes one fewer
% parameter to estimate. 
%=OPTIONAL INPUT
%
%   1) re_bootstrap
%       Recalculate the confidence intervals and probability values. 
function [] = figure4_supp_sum_exp(varargin)
    P = get_parameters;
    if nargin > 0 && varargin{1}
        if exist('Cells', 'var')
            assemble_exp_decay_data(Cells)
        else
            assemble_exp_decay_data()
        end
    else
        load(P.sum_exp_data_path)
    end
    param_list = {'N1', 'alpha', 'kf', 'ks'};
    nparam=numel(param_list);
    figure('Pos', [100, 0, 1200, 1000])
    k = 0;
    n_row = 4;
    n_col = nparam;
    i_label = 1;
    for r = 1:size(T_pval,1)
        idx{1} = strcmp(T_ed.metric, T_pval.metric{r}) & ...
                 strcmp(T_ed.cond_name, T_pval.cond_i{r});
        idx{2} = strcmp(T_ed.metric, T_pval.metric{r}) & ...
                 strcmp(T_ed.cond_name, T_pval.cond_j{r});    
        assert(sum(idx{1})==1)
        assert(sum(idx{2})==1)
        i_label = i_label + 1;
        for p = 1:nparam
            k=k+1;
            ax=subplot(n_row, n_col, k);
            set(ax, P.axes_properties{:}, ...
                    'XLim', [0.5, 2.5], ...
                    'XTick', [])
            for i = 1:2
                x = T_ed.(param_list{p})(idx{i});
                ci = T_ed.([param_list{p} '_CI'])(idx{i}, :);
                err_lower = x - ci(1);
                err_upper = ci(2) - x;
                
                switch T_ed.cond_name(idx{i})
                    case "DV [-10, -2] mm"
                        clr = P.color_order(1,:);
                    case "DV [-2, 0] mm"
                        clr = P.color_order(2,:);
                    case "AP [-8, 0] mm"
                        clr = P.color_order(3,:);
                    case "AP [0, 4] mm"
                        clr = P.color_order(4,:);
                end
                errorbar(i, x, err_lower, err_upper, 'o', 'Color', clr, 'linewidth', 1)
            end
            if r == 1
                title(P.text.(param_list{p}), 'FontWeight', 'Normal')
            end
            if p==1
                label_panel(ax, P.panel_labels(i_label), 'FontSize', P.panel_label_font_size);
            end
            switch param_list{p}
                case {'kf', 'ks'}
                    set(gca, 'YScale', 'log', ...
                             'YLim', [-1e6, -1e-6], ...
                             'YTick', [-1e6, -1, -1e-6], ...
                             'YTickLabel', {'-10^{6}', '-1', '-10^{-6}'})
                otherwise
                    if min(ylim) > 0
                        ylim(ylim.*[0,1]);
                    elseif max(ylim) < 0
                        ylim(ylim.*[1,0]);
                    end
            end
            pval = T_pval.(['pval_2t_', param_list{p}])(r);
            xlabel(['{\it p = ' num2str(pval) '}'], 'FontSize', P.font_size)
        end
    end
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep 'figure4_supp_sum_exp_BtoE'], P.figure_image_format{i})
    end
end