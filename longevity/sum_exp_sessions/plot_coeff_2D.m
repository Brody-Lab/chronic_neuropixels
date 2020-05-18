% PLOT_COEFF_2D plot the coefficient estimtes of the empirical and
% bootstrap samples in 2D. 
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

function [] = plot_coeff_2D()
    %% Parameters
    re_bootstrap = false;
    %% Get files
    P = get_parameters;
    if re_bootstrap
        if exist('Cells', 'var')
            assemble_exp_decay_data(Cells)
        else
            assemble_exp_decay_data()
        end
    end
    load(P.sum_exp_data_path)
    %% Figure
    figure('Pos', [100, 0, 1200, 1440])
    k = 0;
    n_col = 6;
    n_row = 3;
    label_offset = 0;

    param_list = {'alpha', 'ks', 'kf', 'N1'};
    nparam=numel(param_list);
    cond_list = {{'Overall_average'}, ...
                 {'DV [-10, -2] mm', 'DV [-2, 0] mm'}, ...
                 {'AP [-8, 0] mm', 'AP [0, 4] mm'}};
    ncond = numel(cond_list);
    resp_var = ["unit";  "unit/n_elec" ;  "unit/n_elec"];
    colors = {{zeros(1,3)}, ...
              {P.color_order(1,:), P.color_order(2,:)}, ...
              {P.color_order(3,:), P.color_order(4,:)}};
       
    for c=1:ncond
    for i = 1:nparam
        for j = i+1:nparam
            k = k + 1;
            ax=subplot(n_col, n_row,k);
            h = [];
            for cc = 1:numel(cond_list{c})
                h(cc)=plot_error(ax, T_ed, param_list{i}, param_list{j}, colors{c}{cc}, resp_var{c}, cond_list{c}{cc});
            end
            if c < 3;  xlabel(''); end
            if i==1&&j==2
                label_panel(gca, P.panel_labels(c), 'FontSize', P.panel_label_font_size);
                if c == 1
                    text(0.1, 10, '$N=N_{1}[\alpha e^{-(t-1)/\tau_{fast}}) + (1-\alpha)e^{-(t-1)/\tau_{slow}}]$', ...
                         'interpreter', 'latex', 'fontsize', P.font_size*4/5)
                end
%                 legend(h, cond_list{c}, 'location', 'best')
            end
        end
    end
    end
   all_ax = get(gcf, 'children');
%    linkaxes(all_ax([1,4,7]))
%    linkaxes(all_ax([2,5,8]))
%    linkaxes(all_ax([3,6,9]))
   
   %% Save
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep mfilename], P.figure_image_format{i})
    end
end

function [h] = plot_error(ax, T_ed, x, y, clr, resp_var, cond_name)
    idx = strcmp(T_ed.resp_var, resp_var) & ...
          strcmp(T_ed.cond_name, cond_name);
    if sum(idx)~=1
        error('Non-unique condition for %s, %s', resp_var, cond_name)
    end
    xdata = T_ed.([x '_boot']){idx};
    ydata = T_ed.([y '_boot']){idx};
    P = get_parameters;
    set(ax, P.axes_properties{:});
    h=plot(ax, xdata, ydata, 'o', 'Color', clr, 'markersize', 2);
    plot(ax, median(xdata), median(ydata), ...
            '^', 'Color', clr, ...
            'linewidth', 2, ...
            'markersize', 10, ...
            'MarkerfaceColor', 'w');
    switch x
        case {'kf', 'ks'}
            set(gca, 'xscale', 'log', 'xlim',  [-1e6, -1e-6], 'xtick', [-1e6, -1e0, -1e-6])
        case 'alpha'
            set(gca, 'xlim', [0,1])
    end
    switch y
        case {'kf', 'ks'}
            set(gca, 'yscale', 'log', 'ylim', [-1e6, -1e-6], 'ytick', [-1e6, -1e0, -1e-6])
        case 'alpha'
            set(gca, 'ylim', [0,1])
    end
    
    xlabel([P.text.(x)])
    ylabel([P.text.(y)])
end
