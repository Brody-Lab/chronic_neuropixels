% PLOT_GAIN_NOISE_MEDIAN plot the median the median noise level
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = boxplot_noise_stats(varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
addParameter(parseobj,'ylim',[]);
addParameter(parseobj,'var','noise_uV');
addParameter(parseobj,'var2',[]);
addParameter(parseobj,'stat_func',@median);
parse(parseobj, varargin{:});
P_in = parseobj.Results;
analyze_gain_noise_data;
%% plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
q3=norminv(.75);
q=norminv(0.975);
w=(q-q3)/(2*q3); % whisker length for 95% coverage
if ~isempty(P_in.var2)
    tmp=eval(P_in.var);
    if length(tmp(1))==1
        stat_all=eval(P_in.var);
        stat_all2=eval(P_in.var2);         
    else
        stat_all=cellfun(@(x)x(end),eval(P_in.var));
        stat_all2=cellfun(@(x)x(end),eval(P_in.var2));    
    end
    boots_new = bootstrp(1000,P_in.stat_func,cat(1,stat_all{idx_new}),cat(1,stat_all2{idx_new}));
    boots_exp = bootstrp(1000,P_in.stat_func,cat(1,stat_all{~idx_new}),cat(1,stat_all2{~idx_new}));    
else
    stat_all=cellfun(@(x)x(end),eval(P_in.var));
    boots_new = bootstrp(1000,P_in.stat_func,cat(1,stat_all{idx_new}));
    boots_exp = bootstrp(1000,P_in.stat_func,cat(1,stat_all{~idx_new}));
end
box_plot([1 2],[boots_new boots_exp],'center',@mean,'labels',{'New','Expl.'});

if mean(boots_new)>mean(boots_exp)
    big=boots_new;
    small=boots_exp;
else
    big=boots_exp;
    small=boots_new;
end
if prctile(big,2.5)>mean(small) && prctile(small,97.5)<mean(big)
    if prctile(big,0.5)>mean(small) && prctile(small,99.5)<mean(big)
        if prctile(big,0.1)>mean(small) && prctile(small,99.9)<mean(big)
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
    set(gca, P.axes_properties{:}, ...
             'XLim', [0.7, 2.25], ...
             'Xtick', [1 2]);
else
    set(gca, P.axes_properties{:}, ...
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
if strcmp(P_in.var,'noise_uV')
    switch func2str(P_in.stat_func)
        case 'median'
            ylabel({'Median Noise','(\muV_R_M_S)'});
        case '@(x)mean(x>20)*100'
            ylabel('% > 20\muV_R_M_S');        
    end
elseif strcmp(P_in.var,'bank_0_noise_z') && strcmp(func2str(P_in.stat_func),'rsquare')
    ylabel('Across-bank R^2');
end
end