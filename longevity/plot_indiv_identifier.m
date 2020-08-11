function plot_indiv_identifier(T,varargin)
    P=get_parameters;
    p=inputParser;
    p.addParameter('legend_on', true, @(x) isscalar(x) && islogical(x));
    p.addParameter('metric', 'unit', @(x) all(ismember(x, P.longevity_metrics)))
    p.addParameter('normalize_by_electrodes', true, @(x) islogical(x)&&isscalar(x))
    p.parse(varargin{:});
    params=p.Results;
    uidx = unique(T.identifier);
    colors=distinguishable_colors(length(uidx));
    for i=1:length(uidx)
        idx = T.identifier == uidx(i);
        values=T.(params.metric)(idx);
        if params.normalize_by_electrodes
            values = values ./ T.n_elec(idx);
        end
        h(i)=plot(T.days_elapsed(idx),values,'Marker','o','linewidth',1,'color',colors(i,:));hold on;
        [rat(i),serial(i),bank(i)] = decode_identifier(uidx(i));
        legend_string(i) = strcat(rat(i),',',{' '},serial(i),', bank',{' '},bank(i));
    end
    if params.legend_on
        legend(h,legend_string);
    end
    set(gca,P.axes_properties{:},P.custom_axes_properties.longevity{:},'box','off');
    xlabel('Days since implant');
    switch params.metric
        case 'unit'
            ylabel_string='Units';
        case 'single_unit'
            ylabel_string = 'Single units';
        otherwise
            error('not implemented yet.');
    end
    if params.normalize_by_electrodes
        ylabel_string = [ylabel_string,'/electrode'];
    end
    ylabel(ylabel_string);
end
        
    