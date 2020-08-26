%   2016-10-24
%   
%   Set the yaxis of subplots of a figure to have the same limits and tick marks
function [] = regularizeY()

% 1) get information
axes = findobj(gcf, 'Type', 'Axes');
nPlot = length(axes);

yLim = nan(nPlot,2);

for p = 1:nPlot
    yLim(p,:) = get(axes(p), 'YLim');
end

% 2) determine min and max yLimits, as well as tick marks
yLim = [min(yLim(:,1)), max(yLim(:,2))];

% 3) set all plots to have the same yLim
for p = 1:nPlot
    set(axes(p), 'YLim', yLim)
end