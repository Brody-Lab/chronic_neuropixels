% FIGURE6 makes figure 6 from Luo*, Bondy* et al
%
%   The figure show that explanted probes and unimplanted probes have
%   similar input-referred noise and can acquire neural signals of similar
%   quality.
%
%=OPTIONAL INPUT
%
%   from_scratch
%       logical scalar indicating whether the data will be fetched from
%       bucket
function[]=figure6(varargin)
P = get_parameters;
P_input = inputParser;
addParameter(P_input, 'from_scratch', false, @(x) isscalar(x) && islogical(x))
addParameter(P_input,'Cells',[]);
parse(P_input, varargin{:});
P_input = P_input.Results;
if ~isempty(P_input.Cells)
    Cells = P_input.Cells;
end
add_folders_to_path
if P_input.from_scratch
    % note: regnerating these data files from scratch requires raw files and
    % data wrangling code not in the chronic_neuropixels repository.
    get_gain_noise_data
    Cells=collect_cells_files();
    Cells=postprocess_Cells(Cells);  
    add_cumulative_days_implanted
    compute_choice_selectivity    
else
    if ~exist('Cells', 'var')
        fprintf('Loading the variabe CELLS...')
        load([P.data_folder_path filesep 'Cells.mat']);
    end
    if ~exist('ChoiceMod', 'var')    
        fprintf('Loading the variabe ChoiceMod...')        
        load(P.choice_sel_mat_path,'ChoiceMod');
    end
end
figure('Pos', [100, 50, 1500, 800],'color','w')
k = 0;
positions = {[0.05 0.35 0.27 0.62],... %A
             [0.38 0.35 0.27 0.62],... %B
             [0.73 0.76 0.12 0.2],...  %C
             [0.93 0.76 0.06 0.2],...  %D             
             [0.73 0.42 0.12 0.2],...  %E
             [0.93 0.42 0.06 0.2],...  %F             
             [0.05 0.1 0.09 0.2],...   %G
             [0.22 0.1 0.09 0.2],...   %H
             [0.39 0.1 0.09 0.2],...   %I
             [0.56 0.1 0.09 0.2],...   %J
             [0.73 0.1 0.12 0.2],...   %K
             [0.93 0.1 0.06 0.2]};     %L      
         
% A: shank maps of median noise on all explanted probes         
k=k+1;
subplot('position',positions{k});
plot_noise_all_probes('cmap_bounds',[0 50],'stat_func',@(x)min(median(x),50));
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
% add colorbar
c=colorbar('Location','manual',...
    'Position',[0.110 0.8700 0.0250 0.085000],'TickDirection','out','AxisLocation','in');
l=c.Label;
l.String = {'Median','noise','(\muV_R_M_S)'};
l.Rotation=0;
l.Position = [-1.1 128 0];
l.VerticalAlignment='middle';
l.HorizontalAlignment='center';
ticks = [0 25 50]*255/50;
c.Limits=[ticks(1) ticks(end)];
c.Ticks = ticks;
c.TickLabels={'0' '25' '\geq50'};




% B: shank maps of high-noise channels on all explanted probes   
k=k+1;
subplot('position',positions{k});
plot_noise_all_probes('cmap_bounds',[0 .75],'stat_func',@(x)min(0.75,mean(x>20)));
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
% add colorbar
c=colorbar('Location','manual',...
    'Position',[0.4400 0.8700 0.0250 0.085000],'TickDirection','out','AxisLocation','in');
l=c.Label;
l.String = {'Percent >','20\muV_R_M_S'};
l.Rotation=0;
l.Position = [-1.1 128 0];
l.VerticalAlignment='middle';
l.HorizontalAlignment='center';
ticks = [0 .25 .5 0.75]*255/0.75;
c.Limits=[ticks(1) ticks(end)];
c.Ticks = ticks;
c.TickLabels={'0' '25' '50' '\geq75'};


% C: scatter plot of median noise as a function of days implanted
k=k+1;
subplot('position',positions{k});
plot_gain_noise_median('axes', gca)
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% D: box plot of median noise, new v. explanted
k=k+1;
subplot('position',positions{k});
boxplot_noise_stats('ylim',[8 8.4],'var','noise_uV','stat_func',@median,'axes',gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% E: scatter plot of high-noise channels as a function of days implanted
k=k+1;
subplot('position',positions{k});
plot_gain_noise_broken_frac('axes', gca)
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% F: box plot of noisy channels, new v. explanted
k=k+1;
subplot('position',positions{k});
boxplot_noise_stats('ylim',[0 4],'var','noise_uV','stat_func',@(x)mean(x>20)*100,'axes',gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% G/H: plot of number SUs and peak amplitude across time for three
% consecutive implantations of one probe
g=0;
for metric = {'n_good_units','Vpp'}
    k = k + 1;
    g=g+1;
    subplot('position',positions{k});
    plot_recordings_same_probe(Cells, ...
                          'metric', metric{:}, ...
                          'axes', gca, ...
                          'legend', g==1);
   label_hdl(k)= label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
end

% I: choice decoding analysis for the three consecutive implantations
k=k+1;
% plot the average
subplot('position',positions{k});
set(gca, P.axes_properties{:},'xlim',[-1 0.5]);
% avg = cell2mat(cellfun(@(x) mean(abs(x.AUC-0.5)), ChoiceMod, 'uni', 0)');
% hline = plot(ChoiceMod{i}.time_s, avg, 'linewidth', 1);
for i = 1:numel(ChoiceMod)
    hdl(i) = shadedErrorBar(ChoiceMod{i}.time_s, abs(ChoiceMod{i}.AUC-0.5)+0.5, {@mean, @sem});
    hdl(i).mainLine.Color = P.explant_color_order{i};
    hdl(i).mainLine.LineWidth = 1;
    hdl(i).patch.FaceColor = P.explant_color_order{i};
    hdl(i).patch.FaceAlpha = 0.6;
end
yl=ylim;
ylim([0.5 yl(2)]);
plot([0,0],[0.5 yl(2)], 'k-', 'linewidth', 0.5);
xlabel('Time from movement (s)')
ylabel({'Average','choice selectivity'})
label_hdl(k) = label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size)

% J: scatter plots of noise on bank 0 versus bank 1
k=k+1;
subplot('position',positions{k});
bank_noise_scatter('axes', gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% K: scatter plot of across-bank variance explained as a function of days implanted
k=k+1;
subplot('position',positions{k});
plot_noise_bank_corr('axes', gca );
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% L: box plot of across bank R2, new v. explanted
k=k+1;
subplot('position',positions{k});
boxplot_noise_stats('ylim',[0.6 1],'var','bank_0_noise_z','var2','bank_1_noise_z','stat_func',@rsquare,'axes',gca);
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% standardize label position
for i = 1:numel(label_hdl)
    label_hdl(i).Position([2,4]) =  label_hdl(1).Position([2,4]);
    label_hdl(i).Position(1) = positions{i}(1)-0.05;
    if i>6
        label_hdl(i).Position(2)=0.35;    
    elseif i>4
        label_hdl(i).Position(2)=0.7;

    end
end

for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure6'], P.figure_image_format{i})
end