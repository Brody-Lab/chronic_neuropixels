% PLOT_SUM_EXP_TRODES plot the fits to the unit count recorded from
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

function plot_sum_exp_trodes(S, varargin)
    P = get_parameters;
    parseobj = inputParser;
    addParameter(parseobj, 'brain_area', '', ...
        @(x)validateattributes(x,{'char', 'string', 'cell'},{}))
    addParameter(parseobj, 'i_mdl', 1, @(x) validateattributes(x, {'numeric'}, {'integer', 'positive'}));
    parse(parseobj, varargin{:});
    P_in = parseobj.Results;
    yobsv = S.T_trode.(S.P_in.metric);
    t = S.T_trode.days_since_init;
    if isfield(S.SX, 'N1f')
        XN1f = padarray(S.SX.N1f{:,:}, [0,1], 1, 'pre');
        XN1s = padarray(S.SX.N1s{:,:}, [0,1], 1, 'pre');
        ypred = calc_resp_var_N1f_N1s(S.betas, XN1f, XN1s, S.SX.k{:,:}, t);
    else
        XN1 = padarray(S.SX.N1{:,:}, [0,1], 1, 'pre');
        ypred = calc_resp_var_N1_a(S.betas, XN1, S.SX.k{:,:}, t);
    end
    bin_edges = 2.^(-0.5:9.5);
    bin_ctrs = 2.^(0:9);
    
    obsv_avg = groupsummary(yobsv, t+P.x0, bin_edges, ...
                            @nanmean, 'IncludeEmptyGroups', true);
    obsv_sem = groupsummary(yobsv, t+P.x0, bin_edges, ...
                            @sem, 'IncludeEmptyGroups', true);
    pred_avg = groupsummary(ypred, t+P.x0, bin_edges, ...
                            @nanmean, 'IncludeEmptyGroups', true);
    pred_sem = groupsummary(ypred, t+P.x0, bin_edges, ...
                            @sem, 'IncludeEmptyGroups', true);
    figure('Pos', [500,500, 500, 400])
    set(gca, P.axes_properties{:}, ...
             'Xscale', 'log', ...
             'Xtick', bin_ctrs, ...
             'XLim', [min(bin_edges), max(bin_edges)])
    h_obsv=errorbar(bin_ctrs, obsv_avg, obsv_sem, 'k-', 'linewidth', 0.5);
    h_pred=errorbar(bin_ctrs, pred_avg, pred_sem, 'k--', 'linewidth', 0.5);
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