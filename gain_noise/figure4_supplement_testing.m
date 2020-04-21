% figure4_supplement

%=OPTIONAL INPUT
%
%   from_scratch
%       logical scalar indicating whether the data will be fetched from
%       bucket
function[]=figure4_supplement_testing(varargin)
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
figure('Pos', [100, 50, 1500, 800])
k = 0;
n_col = 2;
n_row = 4;
label_x = 0;
label_y = 1.2;

axes_pos_scaling = [1,1,0.9,0.7];

k=k+1;max_z=2;plot_mode='rank';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('Across-bank rank correlation');title(sprintf('abs(z)>%g removed',max_z));
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
k=k+1;max_z=4;plot_mode='rank';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('');title(sprintf('abs(z)>%g removed',max_z));
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
k=k+1;max_z=8;plot_mode='rank';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('');title(sprintf('abs(z)>%g removed',max_z));
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
k=k+1;max_z=Inf;plot_mode='rank';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('');title('no outliers removed');
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
k=k+1;max_z=2;plot_mode='linear';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('Across-bank linear correlation');
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
k=k+1;max_z=4;plot_mode='linear';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('');
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
k=k+1;max_z=8;plot_mode='linear';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('');
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);
k=k+1;max_z=Inf;plot_mode='linear';
subplot(n_col, n_row,k);
plot_noise_bank_corr('axes', gca ,'example_sn',P_input.example_sn,'plot_mode',plot_mode,'max_z',max_z)
set(gca, 'outerpos', get(gca, 'outerpos').*axes_pos_scaling+[0,0.1,0,0],'ylim',[0 1]);ylabel('');
label_hdl(k)=label_panel(gca, P.panel_labels(k), 'FontSize', P.panel_label_font_size);