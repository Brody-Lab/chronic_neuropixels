% PLOT_FIT_SUM_EXP_TRODES plot the fits to the unit count recorded from
% each electrode
%
%=INPUT
%
%   S
%       The structure made by SELECT_MDL_SUM_EXP_TRODES
function plot_fit_sum_exp_trodes(S, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
    addParameter(parseobj, 'i_mdl', 1, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    addParameter(parseobj, 'iteration', 1, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    addParameter(parseobj, 'refit', false, @(x) x==0||x==1);
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    
    i_trodes = S.trodes_used(:, P_in.iteration);
    i_trodes = i_trodes(S.T_trode.days_since_init < P.longevity_time_bin_edges(end));
    
    t = S.T_trode.days_since_init(i_trodes);
    y = S.T_trode.(S.P_in.metric)(i_trodes);
    [XN1, Xk] = partition_T_dsgn(S.T_dsgn, S.T_mdl, P_in.i_mdl, 'i_trodes', i_trodes);
    if P_in.refit
        b = fit_mdl_sum_exp_trodes(XN1, Xk, t, y);
    else
        b = S.T_res.b{P_in.i_mdl}(:, P_in.iteration);
        b = b(~isnan(b));
    end
    
    yhat = predict_y(b, XN1, Xk, t);
    mse = mean((y-yhat).^2);
    fprintf('\nMSE = %0.4f\n', mse)
        
    yobsv = y;
    ypred = yhat;
    bins = discretize(t+P.x0,P.longevity_time_bin_edges);
    for i = 1:max(bins)
        obsv_avg(i) = nanmean(yobsv(bins==i));
        obsv_sem(i) = nanstd(bootstrp(100, @mean, yobsv(bins==i)));
        pred_avg(i) = nanmean(ypred(bins==i));
        pred_sem(i) = nanstd(bootstrp(100, @mean, ypred(bins==i)));
    end
    if isempty(P_in.axes)
        figure('Pos', [500,500, 500, 400])
    else
        axes(P_in.axes)
    end
    set(gca, P.axes_properties{:}, ...
             'Xscale', 'log', ...
             'Xtick', P.longevity_time_bin_centers, ...
             'XLim', [min(P.longevity_time_bin_edges), max(P.longevity_time_bin_edges)])
    h_obsv=errorbar(P.longevity_time_bin_centers, obsv_avg, obsv_sem, 'k-', 'linewidth', 0.5);
    h_pred=errorbar(P.longevity_time_bin_centers, pred_avg, pred_sem, 'k--', 'linewidth', 0.5);
    ylim(ylim.*[0,1])
    
    legend([h_obsv, h_pred], {'Obsv.', 'Pred.'}, 'location', 'northeast')
    xlabel('Days since implant (\it{t})')
    ylabel('Units/electrode (\it{N})')    
    
    unique_t = unique(t);
    for i = 1:numel(unique_t)
        obsv_avg(i) = nanmean(yobsv(unique_t(i)==t));
        obsv_sem(i) = nanstd(bootstrp(100, @mean, yobsv(unique_t(i)==t)));
        pred_avg(i) = nanmean(ypred(unique_t(i)==t));
        pred_sem(i) = nanstd(bootstrp(100, @mean, ypred(unique_t(i)==t)));
        err_avg(i) = mean(abs(yobsv(unique_t(i)==t) - ypred(unique_t(i)==t)));
    end
    figure('Pos', [1e3,500, 500, 400])
    set(gca, P.axes_properties{:}, ...
             'Xscale', 'log', ...
             'Xtick', P.longevity_time_bin_centers)
    h_obsv=errorbar(unique_t+P.x0, obsv_avg, obsv_sem, 'ko-', 'linewidth', 0.5);
    h_pred=errorbar(unique_t+P.x0, pred_avg, pred_sem, 'k*--', 'linewidth', 0.5);
end