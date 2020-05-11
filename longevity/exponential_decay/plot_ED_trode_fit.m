% plot_ED_trode_fit plot the observed and predicted number of units on each
% electode
%
%=INPUT
%
%   S
%       A structure made by SELECT_EXP_MDL
%
function plot_ED_trode_fit(S, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
    addParameter(parseobj, 'i_mdl', 1, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    
%     b = S.T_mdl.b(P_in.i_mdl,:)' .* S.factor_range(:);
    b = S.T_mdl.b(P_in.i_mdl,:)';
    yobsv = S.T_trode.(S.P_in.metric);
    ypred = eval_exp_fit(b, S.T_trode);
    
    t = S.T_trode.days_elapsed + P.x0;
    bins = discretize(t,P.longevity_time_bin_edges);
    for i = 1:max(bins)
        obsv_avg(i) = mean(yobsv(bins==i));
        obsv_sem(i) = std(bootstrp(100, @mean, yobsv(bins==i)));
        pred_avg(i) = mean(ypred(bins==i));
        pred_sem(i) = std(bootstrp(100, @mean, ypred(bins==i)));
    end
    if isempty(P_in.axes)
        figure('Pos', [1e3,500, 500, 400])
    else
        axes(P_in.axes)
    end
    set(gca, P.axes_properties{:}, ...
             'Xscale', 'log', ...
             'Xtick', P.longevity_time_bin_centers, ...
             'XLim', [min(P.longevity_time_bin_edges), max(P.longevity_time_bin_edges)])
    h_obsv=errorbar(P.longevity_time_bin_centers, obsv_avg, obsv_sem, 'k-', 'linewidth', 0.5);
    h_pred=errorbar(P.longevity_time_bin_centers, pred_avg, pred_sem, 'k--', 'linewidth', 0.5);
%     h_obsv=plot(P.longevity_time_bin_centers, obsv_avg, 'ko-', 'linewidth', 0.5);
%     h_pred=plot(P.longevity_time_bin_centers, pred_avg, 'k*--', 'linewidth', 0.5);
    ylim(ylim.*[0,1])
    
    legend([h_obsv, h_pred], {'Obsv.', 'Pred.'}, 'location', 'northeast')
    xlabel('Days since implant (\it{t})')
    ylabel('Units/electrode (\it{N})')    
    
    text(1, 0.5, '$N=N_{1}e^{-(t-1)/\tau}+\varepsilon$', 'interpreter', 'latex', 'fontsize', P.font_size*0.9)
    text(1, 0.3, '$N_{1}=\sum \beta_{AP} AP + ...$', 'interpreter', 'latex', 'fontsize', P.font_size*0.9)
    text(1, 0.15, '$\tau=\sum \beta_{AP} AP + ...$', 'interpreter', 'latex', 'fontsize', P.font_size*0.9)
end

function [yhat] = eval_exp_fit(b, T_trode)
    b=b(:);
    t = T_trode.days_elapsed;
    [X0, Xt] = make_X0_Xt(T_trode, ~isnan(b));
    b = b(~isnan(b));
    yhat = eval_exp(b,X0,Xt,t);
end

function yhat = eval_exp(b,X0,Xt,t)
    n = size(X0,2);
    y0 = X0*b(1:n);
    tDtau = t ./ (Xt*b(n+1:end));
    yhat = y0.*exp(-tDtau);   
end