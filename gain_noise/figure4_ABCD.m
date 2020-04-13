% figure4_ABCD Make panels A-D of figure 4
%
%=OPTIONAL INPUT
%
%   from_scratch
%       logical scalar indicating whether the data will be fetched from
%       bucket
function[]=figure4_ABCD(varargin)
P_input = inputParser;
addParameter(P_input, 'from_scratch', false, @(x) isscalar(x) && islogical(x))
parse(P_input, varargin{:});
P_input = P_input.Results;
add_folders_to_path
if P_input.from_scratch
    get_gain_noise_data
    add_cumulative_days_implanted
end
P = get_parameters;
figure('Pos', [100, 50, 1500, 400])
k = 0;
n_col = 1;
n_row = 4;
label_x = 0;
label_y = 1.2;

axes_pos_scaling = [1,1,0.9,0.7];

k=k+1;
subplot(n_col, n_row,k);
plot_gain_noise_example('unimplanted', 'axes', gca)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0])
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

k=k+1;
subplot(n_col, n_row,k);
plot_gain_noise_example('implanted', 'axes', gca, 'ylabel', false)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0])
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

k=k+1;
subplot(n_col, n_row,k);
plot_gain_noise_median('axes', gca)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling)
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

k=k+1;
subplot(n_col, n_row,k);
plot_gain_noise_broken_frac('axes', gca)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling)
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

% standardize the inner positon
children=get(gcf, 'Children');
pos = cell2mat(arrayfun(@(x) get(x, 'Position'), children, 'uni', 0));
pos(:,2) = mean(pos(:,2)); % lower left corner
pos(:,4)= mean(pos(:,4)); % height
for i = 1:numel(children)
    set(children(i), 'Position', pos(i,:))
end

% standardize label position
for i = 1:numel(label_hdl)
    label_hdl(i).Position([2,4]) =  label_hdl(1).Position([2,4]);
end

for i = 1:numel(P.figure_image_format)
    saveas(gcf, [P.plots_folder_path filesep 'figure4_1'], P.figure_image_format{i})
end