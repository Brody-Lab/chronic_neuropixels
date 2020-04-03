% PLOT_PSYCHOMETRICS plot psychometric curves
%
%=OPTIONAL INPUT
%
%   tethered
%       A scalar specifying whether to plot the untethered (0),
%       untetherd (1)
function [] = plot_psychometrics(varargin)
P_input = inputParser;
addParameter(P_input, 'tethered', 1, @(x) isscalar(x) && (any(x==[0,1])))
parse(P_input, varargin{:});
P_input = P_input.Results;
P = get_parameters;
T = readtable(P.performance_by_rat_path);
T = T(T.tethered == P_input.tethered,:);
x = -30:30;
figure('Position', P.figure_position_psychometrics)
set(gca, P.axes_properties{:}, ...
         'XLim', [min(x), max(x)], ...
         'Ylim', [0, 100]);
plot(xlim, max(ylim)*[1,1], 'k--', 'LineWidth', 0.5)
plot(xlim, 50*[1,1], 'k--', 'LineWidth', 0.5)
plot(0*[1,1], ylim, 'k--', 'LineWidth', 0.5)
for i = 1:size(T,1)
    y = eval_logistic4(x, T.gamma0(i), T.gamma1(i), T.sens(i), T.bias(i));
    plot(x,y, 'LineWidth', 1);
end    
xlabel('#R-#L clicks')
ylabel('Percent chose right')