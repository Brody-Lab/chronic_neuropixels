% figure4_supplement1 - analysis of correlation in RMS noise across banks

%=OPTIONAL INPUT
%
%   from_scratch
%       logical scalar indicating whether the data will be fetched from
%       bucket
function[]=figure4_supplement2(varargin)
P_input = inputParser;
addParameter(P_input, 'from_scratch', false, @(x) isscalar(x) && islogical(x));
addParameter(P_input, 'example_sn',17131311352);
parse(P_input, varargin{:});
P_input = P_input.Results;
add_folders_to_path
if P_input.from_scratch
    get_gain_noise_data
end
P = get_parameters;
figure; set(gcf,'units','normalized','outerposition',[0.38 0.41 0.37 0.55]);
k = 0;
n_col = 1;
n_row = 2;
label_x = 0;
label_y = 1.2;


k=k+1;
subplot(n_col, n_row,k);
plot_noise_all_probes('cmap_bounds',[0 50],'stat_func',@(x)min(median(x),50));
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);

k=k+1;
subplot(n_col, n_row,k);
plot_noise_all_probes('cmap_bounds',[0 1],'stat_func',@(x)mean(x>20));
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
    saveas(gcf, [P.plots_folder_path filesep 'figure4_supplement2'], P.figure_image_format{i})
end