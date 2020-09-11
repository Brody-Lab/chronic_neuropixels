function response_to_minor_comment_13(varargin)
    %% Load data
    from_scratch = false; % Do you want to reassemble the data files from scratch?
    P=get_parameters;
    if from_scratch
        % note: regnerating these data files from scratch requires raw files not in the chronic_neuropixels repository.    
        Cells=collect_cells_files();
        Cells=postprocess_Cells(Cells);
        Cells = add_SNR_data(Cells);
    else
        if ~exist('Cells', 'var')
            fprintf('Loading the variabe CELLS...')
            load([P.data_folder_path filesep 'Cells.mat'])
        end
    end
    
    %%
    recordings_table = read_recordings_log(P.Adrians_recordings_path);    
    agb_files =  (~ismissing(recordings_table.cells_file) & recordings_table.used_in_chronic_npx_ms==1); % find recordings with a cells file path specified and load those files    
    [multiple_probes_recorded,multiple_probes_implanted] = deal(false(size(Cells)));
    event_SNR = cellfun(@(x)x.SNR.event_SNR,Cells);  
    bitscale = cellfun(@(x)x.SNR.bitscale,Cells);    
    bitscale_prm = cellfun(@(x)x.SNR.bitscale_prm,Cells);    
    matching_bitscale = abs(bitscale-bitscale_prm)<0.01;
    rats = recordings_table.rat_name(agb_files);
    multiple_probes_implanted(1:length(rats)) = ismember(rats,{'A230','A241','A243'});
    multiple_probes_recorded(1:length(rats)) = recordings_table.n_probes_recorded(agb_files)>1;  
    h(1)=histogram(event_SNR(~multiple_probes_implanted & matching_bitscale),'NumBins',10,'Normalization','pdf','BinLimits',[0 17],'EdgeColor',[1 1 1]/2,'LineWidth',2,'FaceColor','b','FaceAlpha',0.2);hold on;
    h(2)=histogram(event_SNR(multiple_probes_implanted & matching_bitscale),'NumBins',10,'Normalization','pdf','BinLimits',[0 17],'EdgeColor',[1 1 1]/2,'LineWidth',2,'FaceColor','r','FaceAlpha',0.2);hold on;    
    set(gca,P.axes_properties{:},'box','off');
    legend(h,{sprintf('single-probe implants, n=%d',sum((~multiple_probes_implanted & matching_bitscale))),sprintf('multi-probe implants, n=%d',sum((multiple_probes_implanted & matching_bitscale)))});
    xlabel('Mean Event SNR');
    ylabel('Fraction of Sessions');
    fprintf('Mean event SNR across multi probe sessions is %g.\n',mean(event_SNR(multiple_probes_implanted & matching_bitscale)));
    fprintf('Mean event SNR across single probe sessions is %g.\n',mean(event_SNR(~multiple_probes_implanted & matching_bitscale)));    
    for i = 1:numel(P.figure_image_format)
        saveas(gcf, [P.plots_folder_path filesep 'response_to_minor_comment_13'], P.figure_image_format{i})
    end    
end