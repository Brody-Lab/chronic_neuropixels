% PLOT_AVERAGE_PSYCHOMETRICS plot psychometric curves
%
% OPTIONAL OUTPUT
%
%   1) the axes object
function varargout = plot_average_psychometrics(varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P = get_parameters;
T = readtable(P.performance_by_rat_path);
x = -30:30;
if isempty(P_in.axes)
    figure('Position', P.figure_position_psychometrics)
else
    axes(P_in.axes)
end
set(gca, P.axes_properties{:}, ...
         'XLim', [min(x), max(x)], ...
         'Ylim', [0, 100]);
plot(xlim, max(ylim)*[1,1], 'k--', 'LineWidth', 0.5)
plot(xlim, 50*[1,1], 'k--', 'LineWidth', 0.5)
plot(0*[1,1], ylim, 'k--', 'LineWidth', 0.5)

for tethered = [0,1]
    inds=find(T.tethered == tethered);
    y=nan(numel(inds), numel(x));
    for i = 1:numel(inds)
        ind = inds(i);
        y(i,:) = eval_logistic4(x, T.gamma0(ind), T.gamma1(ind), T.sens(ind), T.bias(ind));        
    end
    wt=sqrt(T.n_trials(inds))/sum(sqrt(T.n_trials(inds)));
    y = y.*wt;
    boots=bootstrp(P.n_boots, @sum, y);
    hdl(tethered+1) = plot(x, mean(boots), 'k-');
    if tethered == 0
        hdl(tethered+1).LineStyle='--';
    end
end
legend(hdl, {'Untethered', 'Tethered'}, 'location', 'best')
xlabel('#R-#L clicks')
ylabel('Percent chose right')
if nargout > 0
    varargout{1} = gca;
end