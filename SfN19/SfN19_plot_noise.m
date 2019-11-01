% 2019_10_12
if ~exist('Res', 'var')
    ans = questdlg('Do you want to re-process the data?');
    if strcmpi(ans, 'yes')
        Res = measure_gain_noise('D:\noise_gain_measurements_2019_10_09', 'max_samples', 1.5e5);
    else
        load('C:\ratter\Analysis\tzluo_data\SfN19_plot_noise\tmp_2019_10_12.mat')
    end
end
%% process the raw results
uV_broken = 30;

S = Res;
% remove probe 17131312042 because I caused a shift register error during
% measurement. *Sigh
S = S([S.probe_sn] ~= 17131312042);

% probe 18194824142 was recorded with calibration on
idx = arrayfun(@(x) x.probe_sn == 18194824142, S);
gain_amp_50 = 1.114688; 
gain_amp_1000 = 1.116029;
S(idx).gain = S(idx).gain/ gain_amp_50;
S(idx).noise_uV_raw = S(idx).noise_uV_raw/gain_amp_1000;
S(idx).noise_uV = S(idx).noise_uV_raw./S(idx).gain;
% Process the results
for i = 1:numel(S)
    S(i).noise_uV([192, 576, 960]) = nan; % reference channels
    S(i).gain([192, 576, 960]) = nan; % reference channels
    S(i).noise_uV(S(i).noise_uV > uV_broken) = uV_broken;
end
chan_by_bank{1} = 1:384;
chan_by_bank{2} = 385:768;
chan_by_bank{3} = 767:960;
kColor = PB_get_constant('Color');

% Import the implant log
Implants = PB_import_implant_log;
for i = 1:numel(S)
    idx = Implants.neuropixels_sn == S(i).probe_sn & ~Implants.acute;
    if sum(idx) ~= 1
        S(i).days_implanted = 0;
        S(i).mm_implanted = 0;
        S(i).channels_implanted = 0;
        warning('Cannot find unique entry in Implants log for probe SN# %i', S(i).probe_sn); 
        continue
    end
    S(i).days_implanted = days(Implants.date_explanted(idx) - Implants.implant_date(idx));
    S(i).mm_implanted = Implants.depth_mm(idx);
    S(i).channels_implanted = floor(S(i).mm_implanted*100); % the depth value in the Implants already had the length of tip subtracted
end
[~, idx ] = sort([S.days_implanted]);
S = S(idx);
% make folder for saving
script_name = mfilename;
plots_fldr = ['C:/ratter/Analysis/tzluo/Plots/' script_name];
if ~exist(plots_fldr, 'dir'), mkdir(plots_fldr), end
%% Plot noise
fig_hdl = wide_figure;
set(fig_hdl, 'Position', [ 0         918        1300         420])
for i = 1:numel(S)
    subplot(1,numel(S), i)
    fig_prepare_axes
    
%     xlabel('Channel')
    yticklabel = cellfun(@num2str, num2cell(0:5:uV_broken), 'uni', 0);
    yticklabel{end} = ['\geq' num2str(uV_broken)];
    set(gca, 'Xlim', [0 961], ...
             'xtick', [], ...
             'YLim', [0, uV_broken], ...
             'Ytick', 0:5:uV_broken, ...
             'YtickLabel', yticklabel)
    if i == 1, ylabel('RMS noise (uV)'), end
    if i > 1, set(gca, 'ytick', []), end
    title(num2str(S(i).probe_sn))
    hdl = [];
    plot(S(i).noise_uV, 'ko', 'MarkerSize', 3);
    for bank = 0:2
        idx = bank+1;
        noise_median = nanmedian(S(i).noise_uV(chan_by_bank{idx}));
        first_last_chans = chan_by_bank{idx}([1,end]);
        hdl(idx) = plot(first_last_chans, noise_median*[1,1], 'Color', kColor.default(idx, :), 'linewidth', 2);
%         area(first_last_chans, uV_broken*[1,1], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', kColor.default(idx, :));
    end
    
    
    hdl_portion_imp = area([0 S(i).channels_implanted], max(ylim)*[1,1], 'FaceAlpha', 0.1, 'EdgeColor', 'none', 'FaceColor', zeros(1,3));
    if S(i).days_implanted == 0
        title('unimplanted')
    else
        title([num2str(S(i).days_implanted) ' days implanted'])
    end
    xlabel('Channels')
    if i == 1
%         legend(hdl, {'bank 0', 'bank 1', 'bank 2'}, 'Location', 'Best', 'box', 'off')
        text(100, 2, 'bank 0', 'FontSize', 16, 'Color', kColor.default(1, :))
        text(584, 2, '1', 'FontSize', 16, 'Color', kColor.default(2, :))
        text(864, 2, '2', 'FontSize', 16, 'Color', kColor.default(3, :))
    elseif i == 2
        legend(hdl_portion_imp, 'Segment implanted', 'Location', 'south', 'fontsize', 12)
    end
end
saveas(gcf, [plots_fldr filesep 'rms_noise_uV.png'])
saveas(gcf, [plots_fldr filesep 'rms_noise_uV.svg'])
%%
return
%% Plot noise

fig_hdl = wide_figure;
for i = 1:numel(S)
    subplot(1,numel(S), i)
    fig_prepare_axes
    for b = 1:3
        histogram(S(i).noise_uV(chan_by_bank{b}), 0:1:30)
    end
    title(num2str(S(i).probe_sn))
    xlabel('RMS noise (uV)')
end
%% Plot gains
fig_hdl = wide_figure;
for i = 1:numel(S)
    subplot(1,numel(S), i)
    plot(S(i).gain, 'k', 'linewidth', 1);
    xlabel('Channel')
    if i == 1, ylabel('gain'), end
    title(num2str(S(i).probe_sn))
end
regularizeY
%% Plot
fig_hdl = figure;
hdl_old = plot(S(1).noise_uV, 'linewidth', 1);
hold on
hdl_new = plot(S(2).noise_uV, 'linewidth', 1);

xlabel('Channel')
set(gca, 'Xlim', [0, 769])
set(gca, 'YLim', [5, 30], ...
         'Ytick', 5:5:30, ...
         'YtickLabel', {'5', '10', '15', '20', '25', '\geq30'})
ylabel('RMS noise (uV)')

legend([hdl_old, hdl_new], {'Recovered', 'Unimplanted'}, 'Location', 'Best')
if ispc
    plots_fldr = ['C:/ratter/Analysis/tzluo/Plots/' mfilename];
else
    plots_fldr = ['\ratter\Analysis\tzluo\Plots\' mfilename];
end
if ~exist(plots_fldr, 'dir'), mkdir(plots_fldr), end
saveas(fig_hdl, [plots_fldr filesep 'practice.png'])