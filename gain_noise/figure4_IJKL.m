% figure4_IJKL Make panels I-L of figure 4
%
% plot the choice selectivity of isolated units recorded from a
% probe implanted in separate animals
%
%=OPTIONAL INPUT
%
%   1) from_scratch
%       A logical scalar specifying whether to make the plots from scratch
function [] = figure4_IJKL(varargin)
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
    plot_choice_mod(ChoiceMod{i}.time_s, ChoiceMod{i}.AUC, ...
                       'ax', gca, ...
                       'color_map_range', [0 1], ...
                       'merge_LR', true, ...
                       'color_bar', false)
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
% avg = cell2mat(cellfun(@(x) mean(abs(x.AUC-0.5)), ChoiceMod, 'uni', 0)');
% hline = plot(ChoiceMod{i}.time_s, avg, 'linewidth', 1);
for i = 1:numel(ChoiceMod)
    hdl(i) = shadedErrorBar(ChoiceMod{i}.time_s, abs(ChoiceMod{i}.AUC-0.5), {@mean, @sem});
    hdl(i).mainLine.Color = P.color_order(i,:);
    hdl(i).mainLine.LineWidth = 1;
    hdl(i).patch.FaceColor = P.color_order(i,:);
    hdl(i).patch.FaceAlpha = 0.3;
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

% show some stats
peak = cellfun(@(x) max(abs(x.AUC-0.5),[],2), ChoiceMod, 'uni', 0)';
days_elapsed=[];
for i = 1:numel(ChoiceMod)
    days_elapsed =[days_elapsed; repmat(i, size(peak{i}))];
end
pval=anovan(cell2mat(peak), days_elapsed, 'display', 'off', 'continuous', 1);
fprintf('\nANOVA on the the peak choice selectivity with the main factor as number of previous implants')
fprintf('\n   p=%f0.3', pval)