% FIGURE4B plot the quality metrics of neural units recorded from a probe
% that was implanted in separate animals
%
%=OPTIONAL INPUT
%
%   1) from_scratch
%       A logical scalar specifying whether to make the plots from scratch
function [] = figure4_2(varargin)
P=get_parameters;
if nargin > 0
    from_scratch = varargin{1};
else
    from_scratch = false;
end
if from_scratch
    collect_cell_files
    postprocess_Cells
else
    if ~exist('Cells', 'var')
        fprintf('Loading the variabe CELLS...')
        load([P.data_folder_path filesep 'Cells.mat'])
    end
end
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