% PLOT_NOISE_ALL_PROBES plots noise for all explanted probes, relative to
% brain surface
%
function [] = plot_noise_all_probes(varargin)
    %% load parameters and set some things up
    P=get_parameters;
    P_in = inputParser;
    addParameter(P_in,'stat_func',@median);
    addParameter(P_in,'axes',gca);
    addParameter(P_in,'cmap_bounds',[0 1]);
    P_in.parse(varargin{:});
    P_in = P_in.Results;
    T= readtable(P.gain_noise_log_path);
    is_explanted = T.date_explanted>datetime('2000-01-01');
    ycoords = cumsum(repmat([20;0],384/2,1));   
    ycoords = [ycoords;ycoords+3840;ycoords(1:192)+3840*2];
    count=0;
    explanted_idx = find(is_explanted);
    [all_noise,all_distances]=deal([]);
    ylim=[-10 6];
    bin_size_um = 100;
    edges = linspace(ylim(1),ylim(2),range(ylim)*1000/bin_size_um +1);
    bin_mids = 0.5*(edges(1:end-1)+edges(2:end));    
    % loop over explantations
    val = NaN((sum(is_explanted)+1)*2,length(edges)-1);
    [~,depth_idx] = sort(T.electrodes_implanted(is_explanted),'descend');    
    depth_idx = [depth_idx(:);depth_idx(end)+1];
    for idx=depth_idx(:)'
        count=count+1;
        if idx>sum(is_explanted)
            %% average across all probes at end   
            bin_idx = discretize(all_distances/1e3,edges);            
            for k=1:(length(edges)-1)
                if any(bin_idx==k)                
                    val(count*2,k) = P_in.stat_func(all_noise(bin_idx==k));        
                end
            end    
        else
            %% load from table
            data_file_path = [P.gain_noise_fldr_path filesep T.recording_id{explanted_idx(idx)} '.csv'];
            D = readtable(data_file_path);    
            distance_from_surface = ycoords - T.electrodes_implanted(explanted_idx(idx))*10;
            all_distances = [all_distances;distance_from_surface(:)];
            all_noise = [all_noise;D.noise_uV(:)];
            bin_idx = discretize(distance_from_surface(:)/1e3,edges);      
            for k=1:(length(edges)-1)
                if any(bin_idx==k)
                    val(count*2,k) = P_in.stat_func(D.noise_uV(bin_idx==k));                                     
                end
            end                 
        end
    end
    %plot
    val = 255*((val - P_in.cmap_bounds(1) + range(P_in.cmap_bounds)*0.01)./(range(P_in.cmap_bounds)*1.01));
    val(isnan(val))=0;
    image((1:(sum(is_explanted)+1)*2),bin_mids,(val'));
    map=viridis(255);
    map(1,:)=[1 1 1];
    colormap(map);
    set(gca,P.axes_properties{:},'YDir','normal','xaxislocation','origin',...
        'xtick',[],'tickdir','out','ygrid','off','ylim',[ylim(1)-0.5 ylim(2)],...
        'ytick',[ylim(1):2:ylim(2)],'yticklabel',abs([ylim(1):2:ylim(2)])); box off    
    text(2,0.7,{'Brain','Surface'},'FontSize',10,'HorizontalAlignment','Center');
    ylabel('            mm below surface              mm above surface');
    count=0;
    exp_count=0;
    exp_string = {'Initial','Second','Third'};
    line([1 1]*sum(is_explanted)*2+0.7,ylim,'linestyle',':','linewidth',1.5,'color',[1 1 1]/2);    
    for i=depth_idx(:)'
        % bank 0/1 border
        count=count+1;
        top=bin_mids(find((val(count*2,:)),1,'last'))+(edges(2)-edges(1))/2;
        bottom=bin_mids(find((val(count*2,:)),1,'first'))-(edges(2)-edges(1))/2;            
        h=patch([-0.5 -0.5 0 0.5 0.5 ]+count*2,[top bottom bottom-0.3 bottom top ],'k');
        h.FaceColor='none';        
        if i<=sum(is_explanted)
            border = +3.84 - T.electrodes_implanted(i)/1e2;
            line([-0.5 0.5]+count*2,[1 1]*border,'color',[1 1 1],'linewidth',2.5,'linestyle',':');
            border = +3.84*2 - T.electrodes_implanted(i)/1e2;
            line([-0.5 0.5]+count*2,[1 1]*border,'color',[1 1 1],'linewidth',2.5,'linestyle',':'); 
            if T.probe_sn(i)==17131312432
                exp_count=exp_count+1;
                text(count*2-0.9,-3,sprintf('%s',exp_string{exp_count}),'Rotation',90,'color',P.explant_color_order{exp_count},'HorizontalAlignment','center');            
            end
        else
            text(count*2-0.9,-6.5,'Combined','Rotation',90,'FontSize',13,'color','k','HorizontalAlignment','center');            
        end
    end
end