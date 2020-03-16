% make figure 2

%% import recordings log
% recordings_log_csv_file = 'C:\Users\abondy\Downloads\recordings_log - All (5).csv';
% recordings_table = read_recordings_log(recordings_log_csv_file);
% has_cells_file = find(~ismissing(recordings_table.cells_file_figure2));
% for i=1:length(has_cells_file)
%     Cells(i) = load(recordings_table.cells_file_figure2(has_cells_file(i)));
% end

n_electrodes = [383 383 191];
for i=1:length(Cells)
    Cells(i) = import_penetration(Cells(i));
    Cells(i).days_since_surgery = days(datetime(Cells(i).sess_date) - datetime(Cells(i).penetration.date_implanted));
    if length(unique(Cells(i).bank))>1
        warning('Multiple banks recorded for rat %s on %s, index %g.',Cells(i).rat,Cells(i).sess_date,i);
    else
        Cells(i).unique_bank = unique(Cells(i).bank);
    end
    out_of_bank = false(length(Cells(i).ks_good),1);
    bank_no = regexprep(Cells(i).mat_file_name,'.*bank([0-9]).*','$1');
    if length(bank_no)~=1
        warning('No bank no in filename, index %g: %s',i,bank_no);
    elseif str2num(bank_no)~=Cells(i).unique_bank
        warning('Bank no in filename (%g) does not match bank no in Cells file (%g): %g, %s',str2num(bank_no),Cells(i).unique_bank,  i,Cells(i).mat_file_name);
    elseif length(unique(Cells(i).bank))>1
        Cells(i).unique_bank = str2num(bank_no);
        out_of_bank = Cells(i).bank~=Cells(i).unique_bank;
    end
    Cells(i).identifier = [Cells(i).rat,Cells(i).probe_serial,num2str(Cells(i).unique_bank)];
    Cells(i).n_units = length(Cells(i).ks_good(~out_of_bank));
    Cells(i).n_good_units = sum(Cells(i).ks_good(~out_of_bank));
    Cells(i).n_electrodes_in_brain = min((Cells(i).penetration.depth_inserted-0.2)*100 - 384*Cells(i).unique_bank,n_electrodes(Cells(i).unique_bank+1));
    Cells(i).mean_DV = mean(Cells(i).DV(Cells(i).ks_good(:) & ~out_of_bank(:)));
end

colors = distinguishable_colors(100);
unique_identifiers = unique({Cells.identifier});
h(1) = figure('Name','Units over time'); % number of total units over time for each unique recording (rat*probe*bank)
h(2) = figure('Name','Good units over time'); % number of single units over time for each unique recording (rat*probe*bank)
h(3) = figure('Name','Yield over time'); % number of single units over time for each unique recording (rat*probe*bank)
h(4) = figure('Name','Depth over time'); % number of single units over time for each unique recording (rat*probe*bank)

clear k j
for i=1:length(unique_identifiers)
    current_sessions{i} = find(strcmp({Cells.identifier},unique_identifiers{i}));
    current_bank = str2num(unique_identifiers{i}(end));
    figure(h(1));
    days_elapsed = [Cells(current_sessions{i}).days_since_surgery]+1;
    [days_elapsed,time_idx] = sort(days_elapsed);
    k(i) = plot(days_elapsed,[Cells(current_sessions{i}(time_idx)).n_units],'Color',colors(i,:),'Marker','o','MarkerFaceColor',colors(i,:),'LineWidth',2);hold on
    if i==length(unique_identifiers)
        legend(k,unique_identifiers);  
        set(gca,'xscale','log','yscale','log');    
        xlabel('Days Since Implant');ylabel('No. SU and MU');                
    end
    figure(h(2));
    j(i) = plot(days_elapsed,[Cells(current_sessions{i}(time_idx)).n_good_units],'Color',colors(i,:),'Marker','o','MarkerFaceColor',colors(i,:),'LineWidth',2);hold on    
    if i==length(unique_identifiers)
        legend(j,unique_identifiers);    
        set(gca,'xscale','log','yscale','log');
        xlabel('Days Since Implant');ylabel('No. SU');        
    end
    figure(h(3));
    l(i) = plot(days_elapsed,[Cells(current_sessions{i}(time_idx)).n_good_units]./[Cells(current_sessions{i}(time_idx)).n_electrodes_in_brain],'Color',colors(i,:),'Marker','o','MarkerFaceColor',colors(i,:),'LineWidth',2);hold on    
    if i==length(unique_identifiers)
        legend(l,unique_identifiers);    
        set(gca,'xscale','log','yscale','log');
        xlabel('Days Since Implant');ylabel('Yield (units per channel)');
    end    
    figure(h(4));
    m(i) = plot(days_elapsed - min(days_elapsed)+1,[Cells(current_sessions{i}(time_idx)).mean_DV] - [Cells(current_sessions{i}(time_idx(1))).mean_DV],'Color',colors(i,:),'Marker','o','MarkerFaceColor',colors(i,:),'LineWidth',2);hold on    
    n_cells{i} = [Cells(current_sessions{i}(time_idx)).n_good_units];
    if i==length(unique_identifiers)
        legend(m,unique_identifiers);    
        set(gca,'xscale','log','yscale','linear');
        xlabel('Days Since First Recording');ylabel('Change in mean Cell Depth (mm below first recording)');
    end        
end

% depth change plot
elapsed = cat(2,m.XData);
delta_depth = cat(2,m.YData);
n_cells = cat(2,n_cells{:});
figure;scatter(log10(elapsed),delta_depth,'SizeData',n_cells);