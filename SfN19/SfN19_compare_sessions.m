analysis_name = mfilename;
clc
clear
rat = 'T176';
fldr = 'C:\ratter\Analysis\tzluo_data\Neuropixels\';
Pharma_log = PB_import_pharmacology_log;
%% Get data
vl_rat = Pharma_log.Rat == rat;
vc_dates = arrayfun(@(x) string(datestr(x, 'yyyy_mm_dd')), Pharma_log.Date(vl_rat));
Cells = {};
Sessions = {};
k = 0;
for i = 1:numel(vc_dates)
    d = dir([fldr rat '_' vc_dates{i} '*']);
    if numel(d) > 1
        k = k + 1;
        [Cells{k,1},Sessions{k,1}] = NP_import_AP(rat, vc_dates{i}, 'mat_file_fldr', fldr);
        Trials{k,1} = PB_make_Trials(Sessions{k});
    end
end
n_sess = numel(Sessions);
pharma = cellfun(@(x) x.pharma.manip, Sessions);

%% Get population decoder
steps_s.clicks_on = -0.5:0.02:1;
steps_s.cpoke_out = -0.5:0.02:0.5;

ref_event = 'clicks_on';
group_by_area = {{'Cg1', 'M2'}, {'PrL'}};

PopDV = {};
for s = 1:n_sess
for g = 1:numel(group_by_area)
    group_by_area{g} = string(group_by_area{g});
    group_by_area{g} = group_by_area{g}(:)';
    idx_cells = any(Cells{s}.cell_area == group_by_area{g},2);
    PopDV{s,g} = PB_compute_PopDV(struct_index(Cells{s}, idx_cells), Trials{s}, 'steps_s', steps_s.(ref_event), 'ref', ref_event, 'latency_s', 0.1);
end
end
%% Plot correlation between the DV and the evidence
kColor = PB_get_constant('color');
kLine = PB_get_constant('line');
figure('Pos', [500          631       378.33          669])
for g = 1:numel(group_by_area)
    subplot(2,1,g)
    fig_prepare_axes
    k = 0;
    hdl = [];
    for manip = {'saline', 'bilateral'}; manip = manip{:};
        vl_manip = pharma == manip;
        steps_s = mean(cell2mat(cellfun(@(x) x.steps_s, PopDV(vl_manip, g), 'uni', 0)));
        rhos = cellfun(@(x) x.corr_DV_clickdiff.avg, PopDV(vl_manip, g), 'uni', 0);
%         for i = 1:numel(rhos)
%             plot(steps_s, rhos{i}, 'Color', kColor.(manip));
%         end
        plot(steps_s, nanmean(cell2mat(rhos)), 'Color', kColor.(manip));
    end
    quickplot('0xy')
    xlabel(['Time from stimulus on (s)'])
    ylabel('corr(\Deltaclicks, DV)')
    title(join(group_by_area{g}))

%     legend(hdl, {'saline', 'bilateral'}, 'location', 'best')
end
regularizeY
%% Plot DV, separately for each quantile of click difference    
n_bins = 4;
figure('Pos', [500       535.67       817.67       764.33])
k = 0;
for g = 1:numel(group_by_area)
for manip = {'saline', 'bilateral'}; manip = manip{:};
    k = k + 1;
    subplot(2,2,k)
    fig_prepare_axes
    set(gca, 'yLim', [-2.5,1.5], ...
             'xLim', [0 1])
    vi_click_diff = [];
    pokedR = [];
    log_odds = [];
    vi_manip = find(pharma == manip);
    vi_manip = vi_manip(:)';
    for i = vi_manip
        vi_click_diff=[vi_click_diff;Trials{i}.click_diff_hz(PopDV{i,g}.vi_trials)];
        pokedR = [pokedR; PopDV{i,g}.pokedR];
        log_odds = [log_odds; PopDV{i,g}.log_odds];
    end          
    cmap = PB_get_my_colors(n_bins*2);
    for choseR = [0,1]
        vi_cd_choice = vi_click_diff(pokedR == choseR);
        edges = quantile(vi_cd_choice, 0:1/n_bins:1);
        edges(1) = -inf;
        edges(end) = inf;
        for b = 1:n_bins
            vi_trials = vi_click_diff >= edges(b) & ...
                        vi_click_diff < edges(b+1) & ...
                        pokedR == choseR;
            hdl = plot(PopDV{1}.steps_s, nanmean(log_odds(vi_trials,:)), 'Color', cmap(choseR*n_bins+b,:));
        end
    end
    quickplot('0xy')
    ylabel('log(P_R/ P_L)')
    xlabel(['Time from stimulus onset (s)'])
    title([char(join(group_by_area{g})) ' ' manip])
end
end
% regularizeY
% saveas(gcf, [analysis_fldr_name filesep 'DV, quantized.png'])
%% Plot a single time slice
t_step = PopDV{1}.steps_s == 0.4;
Color.L = [0,100,200]/255;
Color.R = [200, 0, 0]/255;
figure('Pos', [500       535.67       300       764.33])
k = 0;
manip = 'saline';
vi_manip = find(pharma == manip);
vi_manip = vi_manip(:)';
i = vi_manip(1);
for g = 1:numel(group_by_area)
    clicks_tally = [];
    log_odds = [];
    
    pokedR = [];
    for i = vi_manip
        clicks_tally = [clicks_tally; PopDV{i,g}.clicks_tally.R_minus_L(:,t_step)];
        log_odds = [log_odds; PopDV{i,g}.log_odds(:,t_step)];
        pokedR = [pokedR; PopDV{i,g}.pokedR];
    end  
    k = k + 1;
    subplot(2,1,k)
    fig_prepare_axes
    set(gca, 'Xlim', [-30,30], ....
             'Ylim', [-4,4])
    quickplot('0xy')
    for choice= 'LR'
        idx = pokedR == (choice == 'R');
        plot(clicks_tally(idx), log_odds(idx), 'o', 'Color', Color.(choice), 'MarkerSize', 2);
        rho = num2str(corr(clicks_tally(idx), log_odds(idx), 'rows', 'complete'), '%.2f');
        x_loc = -25 + 40*(choice == 'R');
        text(x_loc, 3.5, ['\rho = ' rho], 'Color', Color.(choice))
    end
end
regularizeY
regularizeX
xlabel('#R - #L clicks')
%% Plot decoding accuracy
figure('Pos', [500          631       378.33          669])
for g = 1:numel(group_by_area)
    subplot(2,1,g)
    fig_prepare_axes
    set(gca, 'YLim', [0.4 0.9])
    plot(xlim, 0.5*[1,1], 'k-', 'Linewidth', 0.5)
    k = 0;
    hdl = [];
    for manip = {'saline', 'bilateral'}; manip = manip{:};
        vl_manip = pharma == manip;
        steps_s = mean(cell2mat(cellfun(@(x) x.steps_s, PopDV(vl_manip, g), 'uni', 0)));
        accuracy = cellfun(@(x) x.Accuracy.avg, PopDV(vl_manip, g), 'uni', 0);
        accuracy = nanmean(cell2mat(accuracy));
%         for i = 1:numel(rhos)
            plot(steps_s, accuracy, 'Color', kColor.(manip));
%         end
    end
    quickplot('0xy')
    xlabel(['Time from ' fix_underscore(ref_event) ' (s)'])
    ylabel('Decoding accuracy')
    title(join(group_by_area{g}))
%     legend(hdl, {'saline', 'bilateral'}, 'location', 'best')
end
regularizeY
%% Check firing rates
t_edges = 0:0.01:1;
t_ctrs = t_edges(1:end-1)+diff(t_edges);
for s = 1:n_sess
    vl_trials = Trials{s}.responded & ~Trials{s}.violated;
    spike_t = cellfun(@(x) x(vl_trials), Cells{s}.spike_time_s.clicks_on, 'uni', 0);
    spike_t = cellfun(@cell2mat, spike_t, 'uni', 0);
    spike_t = cell2mat(spike_t);
    spike_ct{s} = histcounts(spike_t, t_edges) / sum(vl_trials);
end
figure
fig_prepare_axes
for s = 1:n_sess
    plot(t_ctrs, spike_ct{s}, 'Color', kColor.(Sessions{s}.pharma.manip))
end
ylabel('spikes/trial/timebin')