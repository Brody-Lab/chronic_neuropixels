function [] = plot_elab_mdl_coeff()
    P=get_parameters;
    load(P.sum_exp_trodes.elab_mdl_CI_path);
    betas = S.betas_obsv{:,:};
    beta_names = S.betas_obsv.Properties.VariableNames;
    [iA, iB, ik0, ik] = sort_coef_names(beta_names);
    
    figure('Pos', [0, 50, 2560, 600])
    k = 0;
    n_col = 4;
    n_row = 2;

    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_coefficients(S, iA)
    ylabel('Weight fewer $\leftarrow$ units $\rightarrow$ more', 'interpreter', 'latex')    
    
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_coefficients(S, iB)
    regularizeY
    
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_coefficients(S, ik0)
    ylabel('Weight faster $\leftarrow$ decay $\rightarrow$ slower', 'interpreter', 'latex')    
    ax_hdl = gca;
    ax_hdl.OuterPosition(3) = ax_hdl.OuterPosition(3) * 1.2;
    
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_coefficients(S, ik)
   
    %Save
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep 'sum_exp_trodes_elab'], P.figure_image_format{i})
    end
end
%% prep_axes
function [] = plot_coefficients(S, idx)
    P=get_parameters;
    beta_names=S.betas_obsv.Properties.VariableNames(idx);
    xtlabels = cellfun(@lookup_latex_markup_of_parameter, beta_names, 'uni', 0);
    % set up plot
    n = numel(beta_names);
    set(gca, P.axes_properties{:}, ...
             'XLim', [0.5 n+0.5], ...
             'XTick', 1:n, ...
             'XTicklabelRotation', 0, ...
             'XTickLabel', xtlabels, ...
             'TickLabelInterpreter', 'latex')
    plot(xlim, zeros(1,2), 'k-', 'linewidth', 0.5)
    errorbar(1:n, S.betas_obsv{1,idx}, S.lower(idx), S.upper(idx), 'ko')
    ax_hdl = gca;
    ax_hdl.OuterPosition(3) = ax_hdl.OuterPosition(3) * n/7;
    ax_hdl.OuterPosition(4) = ax_hdl.OuterPosition(4)*0.9;
end

%% sort_params
function [iA, iB, ik0, ik] = sort_coef_names(beta_names)
    iA = arrayfun(@(x) contains(x, 'N1f'), beta_names);
    iB = arrayfun(@(x) contains(x, 'N1s'), beta_names);
    ik0 = arrayfun(@(x) contains(x, {'kf', 'ks'}), beta_names);
    ik = arrayfun(@(x) contains(x, 'k_'), beta_names);
    assert(all(iA + iB + ik0 + ik == 1));
end