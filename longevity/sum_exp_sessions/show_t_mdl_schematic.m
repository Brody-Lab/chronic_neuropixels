function [] = show_t_mdl_schematic(Cells)
    P = get_parameters;
    load(P.sum_exp_data_path)
    idx = T_ed.resp_var == "unit" & T_ed.cond_name=="Overall_average";
    assert(sum(idx)==1)
    p = [T_ed.N1(idx), T_ed.alpha(idx), T_ed.kf(idx), T_ed.ks(idx)];
    x = 2.^(0:0.1:13);
    y_hat = sum_2_exp_decay(x,p);
    N1 = T_ed.N1(idx);
    tau_mdl = T_ed.tau_mdl(idx);
    y_hat = y_hat ./N1;
    
    
    T = get_metrics_from_Cells(Cells);
    bootstat = nan(P.longevity_n_boots,max(T.days_bin));
    for j=1:max(T.days_bin)
        idx = T.days_bin==j;
        if sum(idx) > 1
            bootstat(:,j) = bootstrp(P.longevity_n_boots,@mean,T.unit(idx(:)));
        elseif sum(idx)==1
            bootstat(:,j) = metric(idx(:));
        else
            continue
        end
    end
    bootstat = bootstat/N1;
        
    figure('Pos', [500, 500, 1070, 665])
    subplot(2,2,1)
    set(gca, P.axes_properties{:}, 'Xscale', 'log', ...
                                   'XLim', [2^-0.5, max(x)], ...
                                   'XTick', 2.^(0:2:13), ...
                                   'YTick', [0, exp(-1), 1], ...
                                   'YTickLabel', {'0', '1/e', '1'})
    plot(x, y_hat, 'k-', 'linewidth', 1)
    errorbar(P.longevity_time_bin_centers, mean(bootstat), std(bootstat), 'ko')
    
    ylim(ylim.*[0,1]);
    
    plot([min(xlim), tau_mdl], exp(-1)*[1,1], 'k-')
    plot(tau_mdl*[1,1], [0, exp(-1)], 'k-')
    
    xlabel('Days since implant')
    ylabel('Units (normalized)')
    
    label_panel(gca, P.panel_labels(1), 'FontSize', P.panel_label_font_size);
    
    text(125, 0.9, ['\tau_{model} = ' num2str(round(tau_mdl)) ' days'], ...
        'FontSize', P.font_size);
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep 'figure2_supp_sum_exp_A'], P.figure_image_format{i})
    end
end