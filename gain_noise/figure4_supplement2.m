% FIGURE4_SUPPLEMENT2 plots noise for all explanted probes, relative to
% brain surface
%
function [] = figure4_supplement2()
    %% load parameters and set some things up
    P=get_parameters;
    T= readtable(P.gain_noise_log_path);
    is_explanted = T.date_explanted>datetime('2000-01-01');
    ycoords = cumsum(repmat([20;0],384/2,1));   
    ycoords = [ycoords;ycoords+3840;ycoords(1:192)+3840*2];
    figure; set(gcf,'units','normalized','outerposition',[0.38 0.41 0.23 0.58]);
    count=0;
    explanted_idx = find(is_explanted);
    [all_noise,all_distances]=deal([]);
    % loop over explantations
    for idx=1:(sum(is_explanted)+1)
        count=count+1;
        if idx>sum(is_explanted)
            %% average across all probes at end            
            edges = linspace(min(all_distances),max(all_distances),30);
            bin_idx = discretize(all_distances,edges);
            bin_mids = 0.5*(edges(1:end-1)+edges(2:end));
            for i=1:(length(edges)-1)
                vals(:,i) = bootstrp(1000,@(x)min(mean(x),100),all_noise(bin_idx==i));
            end
            in_brain = bin_mids<=0;
            subplot(sum(is_explanted)+1,1,idx);
            clear h;
            for g=1:2
                if g==1
                    h(g)=shadedErrorBar(bin_mids(in_brain)/1e3,vals(:,in_brain),{@mean,@std}); hold on;
                else
                    h(g)=shadedErrorBar(bin_mids(~in_brain)/1e3,vals(:,~in_brain),{@mean,@std});                   
                    h(g).mainLine.Color=[1 1 1]/2;
                end
                h(g).mainLine.LineWidth=1.5;
            end
            text(-9.8,75,sprintf('Average across probes'),'FontSize',11);      
        else
            %% load from table
            i=explanted_idx(idx);
            data_file_path = [P.gain_noise_fldr_path filesep T.recording_id{i} '.csv'];
            D = readtable(data_file_path);    
            distance_from_surface = ycoords - T.electrodes_implanted(i)*10;
            all_distances = [all_distances;distance_from_surface(:)];
            all_noise = [all_noise;D.noise_uV(:)];
            subplot(sum(is_explanted)+1,1,idx);
            in_brain = distance_from_surface<=0;
            %% plot
            for g=1:2
                if g==1
                    plot(distance_from_surface(in_brain)/1e3,min(D.noise_uV(in_brain),100),'.k');hold on   
                else
                    if any(~in_brain)
                        h=plot(distance_from_surface(~in_brain)/1e3,min(D.noise_uV(~in_brain),100),'.');
                        h.MarkerEdgeColor=[1 1 1]/2;
                    end
                end
            end
            text(-9.8,75,sprintf('probe %d\nexplanted %s',T.probe_sn(i),T.date_explanted(i)),'FontSize',8)    
        end
        %% make plots nice
        set(gca,P.axes_properties{:},'yscale','linear','xlim',[-10 6],'ylim',[0 100],'xtick',-10:2:6,'xticklabel',[],...
           'xgrid','on','ytick',[0   100],'yticklabel',{'0'   '\geq100'} );box off
        yl=get(gca,'ylim');line([0 0],yl,'color','k','linewidth',1,'linestyle',':');set(gca,'ylim',yl);
        if idx>sum(is_explanted)
            xlabel('distance from brain surface, mm');
            set(gca,'xticklabel',-10:2:6);
        end
        if i==4
            ylabel({'Input-referred noise','(\muV_R_M_S)'});
        elseif i==1
            text(-5,55,{'In Brain'},'FontSize',16);
            text(2,55,{'Above','Brain'},'FontSize',16,'Color',[0.5 0.5 0.5]);
        end
    end
    %% save figure
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep 'figure4_supplement2'], P.figure_image_format{i});
    end
end