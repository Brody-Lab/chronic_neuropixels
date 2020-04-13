% figure4_EFGH Make panels E-H of figure 4
%
%=OPTIONAL INPUT
%
%   1) from_scratch
%       A logical scalar specifying whether to make the plots from scratch
function [] = figure4_EFGH(Cells)
P=get_parameters;
figure('Pos', [100,100,1500,300])
k = 0;
n_col = 4;
label_offset = 4;
for metric = {'n_units', 'n_good_units', 'fr', 'Vpp'}
    k = k + 1;
    subplot(1, n_col,k);
    plot_recordings_same_probe(Cells, ...
                          'metric', metric{:}, ...
                          'axes', gca, ...
                          'legend', k==1);
    set(gca, 'outerpos', get(gca, 'outerpos').*[1,1,1,0.8]+[0,0.1,0,0])

    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
end
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure4_2'], P.figure_image_format{i})
end