% response to major comment 5, illustrating that manual curation has a minimal effect on yield dynamics
P=get_parameters;
exclude_3A=true;
normalize_by_electrode=true;

%% get the curated data
Cells_curated = collect_cells_files('curated');
Cells_curated=postprocess_Cells(Cells_curated);
T_curated = get_metrics_from_Cells(Cells_curated,'exclude_3A',exclude_3A);

%% get the matching uncurated data
Cells_uncurated = collect_cells_files('uncurated_matching_curated');
Cells_uncurated = postprocess_Cells(Cells_uncurated);
T_uncurated = get_metrics_from_Cells(Cells_uncurated,'exclude_3A',exclude_3A);

%% Set up the figue
figure('Pos', [100, 50, 700, 1250]);
k = 0;
positions = {[0.14 0.8 0.33 0.17],... %A
             [0.64 0.8 0.33 0.17],... %B
             [0.14 0.55 0.33 0.17],...  %C
             [0.64 0.55 0.33 0.17],...  %D             
             [0.14 0.31 0.33 0.17],...  %E
             [0.64 0.31 0.33 0.17],...  %F     
             [0.14 0.07 0.33 0.17],...  %G
             [0.64 0.07 0.33 0.17]};  %H               

%% plot units for uncurated
k=k+1;
subplot('position',positions{k});
plot_indiv_identifier(T_uncurated,'legend_on',false,'metric','unit','normalize_by_electrode',normalize_by_electrode);
text(1,max(ylim)*0.9,'uncurated','FontSize',P.font_size);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);    

%% plot SU, uncurated
k=k+1;
subplot('position',positions{k});
plot_indiv_identifier(T_uncurated,'legend_on',false,'metric','single_unit','normalize_by_electrode',normalize_by_electrode);
text(1,max(ylim)*0.9,'uncurated','FontSize',P.font_size);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);  

%% plot units for curated
k=k+1;
subplot('position',positions{k});
plot_indiv_identifier(T_curated,'legend_on',false,'metric','unit','normalize_by_electrode',normalize_by_electrode);
text(1,max(ylim)*0.9,'curated','FontSize',P.font_size);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);  

%% plot SU, uncurated
k=k+1;
subplot('position',positions{k});
plot_indiv_identifier(T_curated,'legend_on',false,'metric','single_unit','normalize_by_electrode',normalize_by_electrode);
text(1,max(ylim)*0.9,'curated','FontSize',P.font_size);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);  

%% get unit info
uncurated_units = [T_uncurated.unit];
curated_units = [T_curated.unit];
uncurated_single_units = [T_uncurated.single_unit];
curated_single_units = [T_curated.single_unit];
days_elapsed = [T_curated.days_elapsed];
n_elec = [T_uncurated.n_elec];
unit_ratio = curated_units./uncurated_units;
single_unit_ratio = curated_single_units./uncurated_single_units;

%% plot units scatter
k=k+1;
subplot('position',positions{k});
if normalize_by_electrode
    scatter(uncurated_units./n_elec,curated_units./n_elec,40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
    xlabel({'Units/electrode','before curation'});
    ylabel({'Units/electrode','after curation'});    
else
    scatter(uncurated_units,curated_units,40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
    xlabel('Units before curation');
    ylabel('Units after curation');       
end
[r,p] = corrcoef(uncurated_units,curated_units);
stats_label=sprintf('r = %0.2g, p = %0.2g',r(2),p(2));
text(max(xlim)-diff(xlim)*4/5,max(ylim)-diff(ylim)*4/5,stats_label,'FontSize',P.font_size);
set(gca,P.axes_properties{:});h=refline(1,0);h.LineWidth=1;h.Color=[0 0 0];
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);  

%% plot SU scatter
k=k+1;
subplot('position',positions{k});
if normalize_by_electrode
    scatter(uncurated_single_units./n_elec,curated_single_units./n_elec,40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
    xlabel({'Single units/electrode','before curation'});
    ylabel({'Single units/electrode','after curation'});     
else
    scatter(uncurated_single_units,curated_single_units,40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
    xlabel('Single units before curation');
    ylabel('Single units after curation');       
end
[r,p] = corrcoef(uncurated_single_units,curated_single_units);
stats_label=sprintf('r = %0.2g, p = %0.2g',r(2),p(2));
text(max(xlim)-diff(xlim)*4/5,max(ylim)-diff(ylim)*4/5,stats_label,'FontSize',P.font_size);
set(gca,P.axes_properties{:});h=refline(1,0);h.LineWidth=1;h.Color=[0 0 0];
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

%% plot ratio against days implanted, units
k=k+1;
subplot('position',positions{k});

%scatter(([T_uncurated.unit]),[T_curated.unit]./[T_uncurated.unit],40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
%[r,p] = corrcoef(log10([T_uncurated.unit]),[T_curated.unit]./[T_uncurated.unit]);
scatter(days_elapsed,unit_ratio,40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
[r,p] = corrcoef(log(days_elapsed),unit_ratio);
stats_label=sprintf('r = %0.2g, p = %0.2g',r(2),p(2));
text(max(xlim)-diff(xlim)*4/5,max(ylim)-diff(ylim)*4/5,stats_label,'FontSize',P.font_size);
ylabel({'curated/uncurated','unit ratio'});
%xlabel('Units, uncurated');
xlabel('Days since implant');
set(gca,P.axes_properties{:},P.custom_axes_properties.longevity{:},'xscale','log','ylim',ylim.*[0 1]);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

%% plot ratio against days implanted, single units
k=k+1;
subplot('position',positions{k});
%scatter(([T_uncurated.single_unit]),[T_curated.single_unit]./[T_uncurated.single_unit],40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
%[r,p] = corrcoef(log10([T_uncurated.single_unit]),[T_curated.single_unit]./[T_uncurated.single_unit]);
scatter(days_elapsed,single_unit_ratio,40,'markeredgecolor','k','markerfacecolor','none','linewidth',1);
[r,p] = corrcoef(log(days_elapsed),single_unit_ratio);
stats_label=sprintf('r = %0.2g, p = %0.2g',r(2),p(2));
text(max(xlim)-diff(xlim)*4/5,max(ylim)-diff(ylim)*4/5,stats_label,'FontSize',P.font_size);
ylabel({'curated/uncurated','single unit ratio'});
%xlabel('Single units, uncurated');
xlabel('Days since implant');
set(gca,P.axes_properties{:},P.custom_axes_properties.longevity{:},'xscale','log','ylim',ylim.*[0 1]);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% standardize label position
for i = 1:numel(label_hdl)
    label_hdl(i).Position(1) =  label_hdl(i).Position(1)+0.02;
    label_hdl(i).Position(2) = label_hdl(i).Position(2)-0.02;
    if i==6 || i==8
        label_hdl(i).Position(1)=label_hdl(4).Position(1);    
    end
end

for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'response_to_major_comment_5'], P.figure_image_format{i})
end