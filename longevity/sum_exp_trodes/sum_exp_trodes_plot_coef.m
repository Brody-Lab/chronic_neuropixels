% SUM_EXP_TRODES_PLOT_COEF plot coefficients of the regressors in the
% sum-of-exponentials model fit to the unit count from each
% electrode-recording
%
%=INPUT
%   
%   S
%       A structure made using SELECT_EXP_MDL and EST_COEFF_OF_ED_MDL
%
%=OPTIONAL INPUT
%
%   mdl_inds
%       Linear indices of the models to be plotted
%
%   rescale_parameters
%       A scalar logical specifying whether to scale the parameters by the
%       pre-normalized range of the regressors
function [] = sum_exp_trodes_plot_coef(S, varargin)
    P=get_parameters;
    parseobj = inputParser;
    P = get_parameters;
    addParameter(parseobj, 'mdl_inds', 1, ...
        @(x) validateattributes(x, {'numeric'}, { 'integer', 'nonzero'}))
    addParameter(parseobj, 'rescale_parameters', true, @(x) isscalar(x) && islogical(x))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    figure('Pos', [100, 50, 2000, 600])
    k = 0;
    n_col = 4;
    n_row = 2;
    label_offset = 12;
    
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_multiple(ax_hdl, S, 'N1', 'mdl_inds', P_in.mdl_inds, ...
                                   'rescale_parameters', P_in.rescale_parameters);
    ax_hdl.OuterPosition(3) = ax_hdl.OuterPosition(3)*0.82;
    ax_hdl.OuterPosition(4) = ax_hdl.OuterPosition(4)*0.9;
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_multiple(ax_hdl, S, 'alpha', 'mdl_inds',  P_in.mdl_inds, ...
                                      'rescale_parameters', P_in.rescale_parameters);
    ax_hdl.OuterPosition(3) = ax_hdl.OuterPosition(3)*0.22;
    ax_hdl.OuterPosition(4) = ax_hdl.OuterPosition(4)*0.775;
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_multiple(ax_hdl, S, 'k0', 'mdl_inds',  P_in.mdl_inds, ...
                                   'rescale_parameters', P_in.rescale_parameters);
    ax_hdl.OuterPosition(3) = ax_hdl.OuterPosition(3)*0.45;
    ax_hdl.OuterPosition(1) = ax_hdl.OuterPosition(1)+0.1;
    ax_hdl.OuterPosition(4) = ax_hdl.OuterPosition(4)*0.835;
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_multiple(ax_hdl, S, 'k', 'mdl_inds',  P_in.mdl_inds, ...
                                  'rescale_parameters', P_in.rescale_parameters);
    ax_hdl.OuterPosition(3) = ax_hdl.OuterPosition(3)*0.9;
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    LL = mean(-S.T_res.MSE,2);
    LL_rng = min(LL) - max(LL);
    set(ax_hdl, P.axes_properties{:})
    set(gca, 'xscale', 'log', ...
             'xlim', [0.75, size(S.T_res,1)], ...
             'xtick', [1, 5, size(S.T_res,1)], ...
             'ylim', [LL_rng*max(S.T_res.err_norm), 0], ...
             'ticklabelinterpreter', 'latex')
%     set(gca, 'ytick', [2.32, 2.44])
    y = S.T_res.err_norm * LL_rng;
    plot(6:numel(y), y(6:end),'o', 'Color', 0.6*[1,1,1], 'linewidth', 0.5)
    plot(1, y(1),'o', 'Color', [0.15, 0.15, 0.8], 'linewidth', 0.5, 'markerfacecolor', [0.15, 0.15, 0.8])
    plot(2:5, y(2:5),'ko', 'linewidth', 1)
    xlabel('Regressor subset', 'interpreter', 'latex')
    ylabel('Rel. LL', 'interpreter', 'latex')
    ax_hdl.OuterPosition(4) = ax_hdl.OuterPosition(4)*0.4;
    ax_hdl.OuterPosition(3) = ax_hdl.OuterPosition(3)*0.6;
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    set(ax_hdl, P.axes_properties{:}, 'xtick', [], 'ytick', [])
    text(0.1,0.9, '$N=N_{1}[\alpha e^{k_{fast}(t-1)}+(1-\alpha)e^{k_{slow}(t-1)}]e^{k}$', 'Fontsize', 14, 'interpreter', 'latex')
    text(0.1,0.7, '$N_{1}=\beta_{0}^{N_{1}}+\beta_{AP}^{N_{1}}AP+$...', 'Fontsize', 14, 'interpreter', 'latex')
    text(0.1,0.5, '$k=\beta_{AP}^{k}AP+$...', 'Fontsize', 14, 'interpreter', 'latex')
    %% Save
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep 'figure2_mdl'], P.figure_image_format{i})
    end
end
%% plot_multiple
% Plot multiple model's coefficients simultaneously
%
%=INPUT
%
%   ax
%       The axes object in which the plot is made
%
%   S
%       The structure with all the model data
%
%   param_type
%       The type of parameter to be ploted

%
%=OPTIONAL INPUT
%
%   mdl_inds
%       The model indices to be plotted
%
%   rescale_parameters
%       A scalar logical specifying whether to scale the parameters by the
%       pre-normalized range of the regressors

function [] = plot_multiple(ax, S, param_type, varargin)
    P = get_parameters;
    assert(isa(ax, 'matlab.graphics.axis.Axes'))
    assert(isscalar(ax))
    assert(any(ismember({'alpha', 'k0', 'N1', 'k'}, param_type)));
    parseobj = inputParser;
    P = get_parameters;
    addParameter(parseobj, 'mdl_inds', 1, ...
        @(x) validateattributes(x, {'numeric'}, { 'integer', 'nonzero'}))
    addParameter(parseobj, 'rescale_parameters', true, @(x) isscalar(x) && islogical(x))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    if P_in.rescale_parameters
        [b, cil, ciu] = rescale_parameters(S.T_res, S.T_mdl, S.T_regressor);
    else
        b = S.T_res.b_med;
        cil=S.T_res.cil;
        ciu=S.T_res.ciu;
    end
    b(isnan(b))=0; % parameters that weren't fitted were set as nan
    cil(isnan(cil))=0;
    ciu(isnan(ciu))=0;
    parameter_types = ["alpha", "k0", "k0", ...
        string(cellfun(@get_parameter_type, ...
                        S.T_mdl.Properties.VariableNames, 'uni', 0))];
    b = b(:, parameter_types == param_type);
    cil = cil(:, parameter_types == param_type);
    ciu = ciu(:, parameter_types == param_type);
    parameter_names = ["alpha", "k_fast", "k_slow",...
                       string(S.T_mdl.Properties.VariableNames)];
    parameter_names = parameter_names(parameter_types == param_type);
    xtlabels = cellfun(@lookup_latex_markup_of_parameter, parameter_names, 'uni', 0);
    % set up plot
    n = size(b,2);
    set(ax, P.axes_properties{:}, ...
             'XLim', [0.5 n+0.5], ...
             'XTick', 1:n, ...
             'XTicklabelRotation', 0, ...
             'XTickLabel', xtlabels, ...
             'TickLabelInterpreter', 'latex')
    h=refline(0,0); 
    set(h, 'Color', 'k')
    % plot each model
    for i = 1:numel(P_in.mdl_inds)
        ind=P_in.mdl_inds(i);
        if i == 1
            col = [0.15, 0.15, 0.8];
            lw = 1;
        else
            col = 'k';
            lw = 0.5;
        end
        x=(1:n);
        if numel(P_in.mdl_inds)>1
            x=x-0.3 + i/numel(P_in.mdl_inds)*0.6;
        end
        hdl(i) = errorbar(x, b(ind,:), b(ind,:)-cil(ind,:),  ciu(ind,:)-b(ind,:), 'o', 'Color', col, 'Markersize', 3, 'linewidth',lw);
    end
    if any(P_in.mdl_inds==1) && strcmp(param_type, 'k')
        legend(hdl(1), 'Highest LL', 'Location', 'Best', 'interpreter', 'latex')
    end
    
    ylabeltxt.alpha = '';
    ylabeltxt.k = ['Weight ($day^{-1}$)' newline 'faster $\leftarrow$ decay $\rightarrow$ slower'];
    ylabeltxt.N1 = ['Weight ' newline 'fewer $\leftarrow$ units $\rightarrow$ more'];
    ylabeltxt.k0 = ['Change rate ($day^{-1}$)' newline 'faster $\leftarrow$ decay $\rightarrow$ slower'];
    ylabel(ylabeltxt.(param_type), 'interpreter', 'latex')
    
    switch param_type
        case 'alpha'
            title(['Fraction' newline 'with the' newline 'faster' newline 'decay rate'], 'fontweight', 'normal', 'interpreter', 'latex')
        case 'k0'
            title(['Baseline' newline 'change' newline 'rates'], 'fontweight', 'normal', 'interpreter', 'latex')
        case 'N1'
            title(['Coefficients in the' newline 'initial unit count ($N_{1}$)'], 'fontweight', 'normal', 'interpreter', 'latex')
%             ylim([-0.5, 1.5])
        case 'k'
            title('Coefficients in the change rate ($k$)', 'fontweight', 'normal', 'interpreter', 'latex')
%             ylim([-0.025, 0.015])
            ytl = arrayfun(@num2str, yticks, 'uni', 0);
            yticklabels(ytl)
        case 'k0'
%             ylim([-1.5, 0.5])
    end
end