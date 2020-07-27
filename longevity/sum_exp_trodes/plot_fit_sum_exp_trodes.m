% PLOT_FIT_SUM_EXP_TRODES plot the fits to the unit count recorded from
% each electrode
%
%=INPUT
%
%   S
%       The structure made by SELECT_MDL_SUM_EXP_TRODES
%
%=OPTIONAL INPUT
%
%   brain_area
%       A char vector, string, or cell array specifying the brain areas
%       that are included. 
%
%   i_mdl
%       The i-th best model to be examined.
%
%   iteration
%       The model is fitted across multiple repetitions, each time using a
%       random half of the trials. This specifies which iteratation to use.
%
%   refit
%       A logical scalar to refit the model

function plot_fit_sum_exp_trodes(S, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'brain_area', '', ...
        @(x)validateattributes(x,{'char', 'string', 'cell'},{}))
    addParameter(parseobj, 'i_mdl', 1, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    addParameter(parseobj, 'iteration', 1, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    addParameter(parseobj, 'refit', false, @(x) x==0||x==1);
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    
    i_trodes = S.trodes_used(:, P_in.iteration);
    i_trodes = i_trodes & S.T_trode.days_since_init < P.longevity_time_bin_edges(end);
    if ~isempty(P_in.brain_area)
        i_trodes = i_trodes & ismember(S.T_trode.brain_area, string(P_in.brain_area));
    end
    t = S.T_trode.days_since_init(i_trodes);
    yobsv = S.T_trode.(S.P_in.metric)(i_trodes);
    [XN1, Xk] = partition_T_dsgn(S.T_dsgn, S.T_mdl, P_in.i_mdl, 'i_trodes', i_trodes);
    if P_in.refit
        b = fit_mdl_sum_exp_trodes(XN1, Xk, t, yobsv);
    else
        b = S.T_res.b{P_in.i_mdl}(:, P_in.iteration);
        b = b(~isnan(b));
    end
    
    ypred = predict_y(b, XN1, Xk, t);
    mse = mean((yobsv-ypred).^2);
    fprintf('\nMSE = %0.4f\n', mse)
    obsv_avg = groupsummary(yobsv, t+P.x0, P.longevity_time_bin_edges, ...
                            @nanmean, 'IncludeEmptyGroups', true);
    obsv_sem = groupsummary(yobsv, t+P.x0, P.longevity_time_bin_edges, ...
                            @sem, 'IncludeEmptyGroups', true);
    pred_avg = groupsummary(ypred, t+P.x0, P.longevity_time_bin_edges, ...
                            @nanmean, 'IncludeEmptyGroups', true);
    pred_sem = groupsummary(ypred, t+P.x0, P.longevity_time_bin_edges, ...
                            @sem, 'IncludeEmptyGroups', true);
    figure('Pos', [500,500, 500, 400])
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
    
    %%
    [G,unique_t] = findgroups(t);
    obsv_avg = groupsummary(yobsv, G, @nanmean, 'IncludeEmptyGroups', true);
    obsv_sem = groupsummary(yobsv, G, @sem, 'IncludeEmptyGroups', true);
    pred_avg = groupsummary(ypred, G, @nanmean, 'IncludeEmptyGroups', true);
    pred_sem = groupsummary(ypred, G, @sem, 'IncludeEmptyGroups', true);
    figure('Pos', [1e3,500, 500, 400])
    set(gca, P.axes_properties{:}, ...
             'Xscale', 'log', ...
             'Xtick', P.longevity_time_bin_centers)
    h_obsv=errorbar(unique_t+P.x0, obsv_avg, obsv_sem, 'ko-', 'linewidth', 0.5);
    h_pred=errorbar(unique_t+P.x0, pred_avg, pred_sem, 'k*--', 'linewidth', 0.5);
end