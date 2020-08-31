% FIGURE2_SUPP1 Make figure 4--supplement 4 from the manuscrip written by
% Luo, Bondy, et al.
%
%   The figure show that no degradation of spiking signals was detected
%   over two months in rat medial prefrontal cortex (mPFC), the brain
%   region in which the stability of spiking signals was examined in (Jun
%   et al., 2017)
%% Load Cells
from_scratch = false; % Do you want to reassemble the data files from scratch?
P=get_parameters;
if from_scratch
    Cells=collect_cells_files();
    Cells=postprocess_Cells(Cells);
else
    if ~exist('Cells', 'var')
        fprintf('Loading the variabe CELLS...')
        load([P.data_folder_path filesep 'Cells.mat'])
    end
end
%% Make the figue
figure('Pos', [100, 50, 2500, 600])
k = 0;
n_col = 5;
n_row = 2;
label_offset = 0;
T = get_metrics_from_Cells(Cells, 'condition_on', 'brain_area', ...
                                  'brain_area', {{'other'},{'PrL', 'MO'}});
for metric = {'unit', 'single_unit', 'event_rate', 'Vpp'}; metric=metric{:};
    k = k + 1;
    subplot(n_row, n_col,k);
    hdl = plot_average_stability(T, 'metric', metric, ...
                                   'print_sample_size', mod(k,n_col)==1, ...
                              'axes', gca, ...
                              'legend_on', mod(k,n_col)==1);
    hdl(1).mainLine.Color=zeros(1,3);
    hdl(1).patch.FaceColor=zeros(1,3);
    label_panel(gca, P.panel_labels(k+label_offset), 'FontSize', P.panel_label_font_size);
    title(P.text.(metric))
    lgd=findobj(gcf, 'Type', 'Legend');
    lgd.String = {'mPFC', 'Other regions'};
end
% plot histogram
k = k + 1;
subplot(n_row, n_col,k)
make_Vpp_histogram(T, 'color_order', [0,0,0; P.color_order(2,:)])
%% Save
for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep mfilename], P.figure_image_format{i})
end