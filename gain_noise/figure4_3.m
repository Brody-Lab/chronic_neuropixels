% FIGURE 4C plot the choice selectivity of isolated units recorded from a
% probe implanted in separate animals
%
%=OPTIONAL INPUT
%
%   1) from_scratch
%       A logical scalar specifying whether to make the plots from scratch
function [] = figure4c(varargin)
P = get_parameters;
if nargin > 0
    from_scratch = varargin{1};
else
    from_scratch = false;
end
if from_scratch
    compute_choice_selectivity
end
load(P.choice_sel_mat_path,'ChoiceMod')
figure('Pos', [100,100,1500,350])
k = 0;
title_texts = {'Initial implant';
               'Second implant';
               'Third implant'};
label_offset = 8;
n_col = numel(P.choice_sel_file_names)+1;
for i = 1:numel(P.choice_sel_file_names)
    k = k + 1;
    subplot(1,n_col, k)
    set(gca, P.axes_properties{:})
    PB_plot_choice_mod(ChoiceMod{i}.time_s, ChoiceMod{i}.AUC, ...
                       'ax', gca, ...
                       'color_map_range', [0 1], ...
                       'merge_LR', true, ...
                       'color_bar', i==1)
    plot([0,0],ylim, 'w-', 'linewidth', 1)
        set(gca, 'outerpos', get(gca, 'outerpos').*[1,0,1,0.8]+[0,0.05,0,0])
    xlabel('Time from movement (s)')
    title(title_texts{i}, 'fontweight', 'normal', 'Color', P.color_order(i,:))
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size)
end

% plot the average
k = k + 1;
subplot(1,n_col, k)
set(gca, P.axes_properties{:})
avg = cell2mat(cellfun(@(x) mean(abs(x.AUC-0.5)), ChoiceMod, 'uni', 0)');
hline = plot(ChoiceMod{i}.time_s, avg);
for i = 1:numel(hline)
    hline(i).Color = P.color_order(i,:);
end
ylim(ylim.*[0,1]);
plot([0,0],ylim, 'k-', 'linewidth', 0.5)
xlabel('Time from movement (s)')
ylabel('Average choice selectivity')
set(gca, 'outerpos', get(gca, 'outerpos').*[1,0,1,0.8]+[0.02,0.05,0,0])
label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size)

%standardize sizes
children=get(gcf, 'Children');
pos = cell2mat(arrayfun(@(x) get(x, 'Position'), children, 'uni', 0));
pos(:,2) = mean(pos(:,2)); % lower left corner
pos(:,4)= mean(pos(:,4)); % height
for i = 1:numel(children)
    set(children(i), 'Position', pos(i,:))
end


for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure4_3'], P.figure_image_format{i})
end
