%=INPUT
%   
%   S
%       A structure made using SELECT_EXP_MDL and EST_COEFF_OF_ED_MDL
%
%=OPTIONAL INPUT
%
%   mdl_inds
%       Linear indices of the models to be plotted in panels N and O
function [] = plot_ED_trode_results(S, varargin)
    P=get_parameters;
    parseobj = inputParser;
    P = get_parameters;
    addParameter(parseobj, 'mdl_inds', 1:10, ...
        @(x) validateattributes(x, {'numeric'}, { 'integer', 'nonzero'}))
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    figure('Pos', [100, 50, 1300, 600])
    k = 0;
    n_col = 2;
    n_row = 2;
    label_offset = 12;
    for make_this_plot = []
        % plot the regressor coefficients for factors whose units are day^-1 *
        % mm^-1 for the initial value
        k = k + 1;
        ax_hdl=subplot(n_row,n_col, k);
        make_plot(ax_hdl, S, 'y0');
        ylabel('Weight')
        title('Coefficients in the initial value', 'fontweight', 'normal')
        label_panel(ax_hdl, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
        ylim([-0.6,1])
        % plot the regressor coefficients for factors whose units are day^-1 for
        % the initial value
        k = k + 1;
        ax_hdl=subplot(n_row,n_col, k);
        make_plot(ax_hdl, S, 'k');
        ylabel('Weight (day)')
        title('Coefficients in the time constant', 'fontweight', 'normal')
        label_panel(ax_hdl, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
        % plot the frequency some regressors show up in the top 25 models
    
        n = numel(P.ED_trode_regressors_symbols);
        varnames = P.ED_trode_regressors_symbols([2:n/2, n/2+2:end]);
        n = numel(varnames);
        varnames = cellfun(@(x) ['$' x '$'], varnames, 'uni', 0);
        k = k + 1;
        ax_hdl=subplot(n_row,n_col, k);
        set(ax_hdl, P.axes_properties{:}, ...
                 'XLim', [0.5, n+0.5], ...
                 'XTick', 1:n, ...
                   'XTicklabelRotation', 0, ...
                   'TickLabelInterpreter', 'latex', ...
                 'XTickLabel', varnames)
        freq=sum(S.T_mdl{1:25, P.ED_trode_regressors})/25;
        hdl = bar(1:2:n, freq(1:2:end), 0.4, 'k');
        set(hdl, 'faceColor', 0.6*[1,1,1])
        hdl = bar(2:2:n, freq(2:2:end), 0.4, 'k');
        set(hdl, 'faceColor', 0.9*[1,1,1])
        title('Frequency in the 25 models with lowest MSE', 'fontweight', 'normal')
        label_panel(ax_hdl, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    end
    %%
%     k = k + 1;
%     ax_hdl=subplot(n_row,n_col, k);
%     plot_ED_trode_fit(S, 'axes', ax_hdl)
%     label_panel(ax_hdl, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
%     ax_hdl.Position(3)=ax_hdl.Position(3)*0.75;
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_multiple(ax_hdl, S, 'y0', P_in.mdl_inds);
    label_panel(ax_hdl, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    %%
    k = k + 1;
    ax_hdl=subplot(n_row,n_col, k);
    plot_multiple(ax_hdl, S, 'k', P_in.mdl_inds);
    label_panel(ax_hdl, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    %% Save
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep 'figure2_mdl'], P.figure_image_format{i})
    end
end
function [] = make_plot(ax_hdl, S, param_type)
    P = get_parameters;
    idx = contains(P.ED_trode_regressors_all, param_type);
    n = sum(idx);
    x=1:n;
    set(ax_hdl, P.axes_properties{:}, ...
             'XLim', [0.5 n+0.5], ...
             'XTick', x, ...
             'XTicklabelRotation', 0, ...
             'XTickLabel', P.ED_trode_regressors_symbols(idx), ...
             'TickLabelInterpreter', 'latex')
    h=refline(0,0); 
    set(h, 'Color', 'k')
    
    b = S.T_mdl.b(2,idx);
    cil = S.T_mdl.cil(2,idx);
    ciu = S.T_mdl.ciu(2,idx);
        
    if strcmp(param_type, 'k')
        b = asinh(b);
        cil=asinh(cil);
        ciu=asinh(ciu);
        odd_color = 0.9*[1,1,1];
        even_color = 0.6*[1,1,1];
    else
        odd_color = 0.6*[1,1,1];
        even_color = 0.9*[1,1,1];
    end
    
    odd =1:2:n;
    even = 2:2:n;
    hdl = bar(odd, b(odd),  0.4, 'k');
    set(hdl, 'faceColor', 0.6*[1,1,1])
    hdl = bar(even, b(even),  0.4, 'k');
    set(hdl, 'faceColor', 0.9*[1,1,1])   
    
    errorbar(x, b, b-cil,  ciu-b, 'ko',  'markersize', 0.1, 'linewidth',1)
    
    if strcmp(param_type, 'k')
        yt = [-1e4, -1e2, 0, 1e2, 1e4];
        asinh_yt = asinh(yt);
        yticks(asinh_yt);
        ylim([min(asinh_yt), max(asinh_yt)])
        ytl = cellfun(@num2str, num2cell(yt), 'uni', 0);
        yticklabels(ytl)
    end
end
function [] = plot_multiple(ax, S, param_type, inds)
    P = get_parameters;
    idx = contains(P.ED_trode_regressors_all, param_type);
    n = sum(idx);
    set(ax, P.axes_properties{:}, ...
             'XLim', [0.5 n+0.5], ...
             'XTick', 1:n, ...
             'XTicklabelRotation', 0, ...
             'XTickLabel', P.ED_trode_regressors_symbols(idx), ...
             'TickLabelInterpreter', 'latex')
    h=refline(0,0); 
    set(h, 'Color', 'k')
    
    for i = 1:numel(inds)
        ind=inds(i);
        b = S.T_mdl.b(ind,idx);
        b(isnan(b))=0;
        cil = S.T_mdl.cil(ind,idx);
        ciu = S.T_mdl.ciu(ind,idx);
        
        if strcmp(param_type, 'k')
            b = asinh(b);
            cil=asinh(cil);
            ciu=asinh(ciu);
            odd_color = 0.9*[1,1,1];
            even_color = 0.6*[1,1,1];
        else
            odd_color = 0.6*[1,1,1];
            even_color = 0.9*[1,1,1];
        end
        
        x=(1:n)-0.3 + i/numel(inds)*0.6;
        errorbar(x, b, b-cil,  ciu-b, 'ko',  'markersize', 3, 'linewidth',0.5)
    end
    if strcmp(param_type, 'k')
        yt = [-1e4, -1e3, -1e2, -1e1 0, 1e1 1e2, 1e3, 1e4];
        
        asinh_yt = asinh(yt);
        yticks(asinh_yt);
        ylim([min(asinh_yt), max(asinh_yt)])
        ytl = cellfun(@num2str, num2cell(yt), 'uni', 0);
        ytl=arrayfun(@(x) sprintf('$10^%i$',x), log10(abs(yt)), 'uni', 0);
        ytl(yt==0) = {'0'};
        ytl(yt<0) = cellfun(@(x) ['-' x ''], ytl(yt<0), 'uni', 0);
        yticklabels(ytl')
        title('Coefficients in the time constant', 'fontweight', 'normal')
        ylabel('Weight (day)')

    else
        title('Coefficients in the initial value', 'fontweight', 'normal')
        ylabel('Weight')
        ylim([-0.6,1])
    end
end