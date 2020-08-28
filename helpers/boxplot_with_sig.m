% PLOT_GAIN_NOISE_MEDIAN plot the median the median noise level
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = boxplot_with_sig(var,var2,varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', gca, @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj,'ylim',[]);
addParameter(parseobj,'stat_func',@mean);
addParameter(parseobj,'labels',{'' ''});
addParameter(parseobj,'ylabel','');
parse(parseobj, varargin{:});
P_in = parseobj.Results;
%% plot
axis(P_in.axes);
q3=norminv(.75);
q=norminv(0.975);
w=(q-q3)/(2*q3); % whisker length for 95% coverage
n_boots=1e4;
boots1 = bootstrp(n_boots,P_in.stat_func,var);
boots2 = bootstrp(n_boots,P_in.stat_func,var2);
box_plot([1 2],[boots1 boots2],'center',@mean,'labels',P_in.labels);
if mean(boots1)>mean(boots2)
    big=boots1;
    small=boots2;
else
    big=boots2;
    small=boots1;
end
p_val = (sum(small>big)+1)/n_boots;
if p_val<0.05
    if p_val<0.01
        if p_val<0.001
            %p<0.001
            n_stars=3;
        else%p<0.01
            n_stars=2;
        end
    else % p<0.05
        n_stars=1;
    end
else
    %not sig
    n_stars=0;
end
  


if isempty(P_in.ylim)
    set(gca, ...
             'XLim', [0.7, 2.25], ...
             'Xtick', [1 2]);
else
    set(gca,  ...
                 'XLim', [0.7, 2.25], ...
                 'Xtick', [1 2], ...
                 'ylim',P_in.ylim);    
end
box off;
if n_stars
    yl=get(gca,'ylim');
    bottom = prctile(big,97.5)+range(get(gca,'ylim'))/10;
    height = range(get(gca,'ylim'))/20;
    if height+bottom>yl(2)
        bottom = 0.5*plus(yl(2),prctile(big,97.5));        
        height = yl(2) - bottom;
    end
    for i=1:2
        line([i i],[0 height]+bottom,'linewidth',1,'color','k');
    end
    line([1 2],[1 1]*height + bottom,'linewidth',1,'color','k');
    text(1.5,bottom+height*1.5,repmat('*',1,n_stars),'FontSize',15,'HorizontalAlignment','center');
end
ylabel(P_in.ylabel);
end