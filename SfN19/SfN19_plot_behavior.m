% SfN19_plot_behavior
%
% Show the psychometric of animals recorded from during during behavior in
% the Poisson Clicks task
clc
clear
script_name = mfilename;
%% Specify rats and sessions, tethered only
% I am excluding T170, T182, T196, and T209 because their behavioral
% performance degraded severely after the surgery

S = struct;

% T176
Pharma = PB_import_pharmacology_log;
idx = Pharma.Rat == 'T176' & Pharma.Manipulation == 'saline';
S.rat = repmat("T176", sum(idx),1);
S.date = Pharma.Date(idx);
S.tethered = true(sum(idx),1);

% T173
n = numel(S.date);
S.date(n+1) = datetime('2018-04-16');
S.rat(n+1) = 'T173';
S.tethered(n+1) = true;

%K265
n = numel(S.date);
S.date(n+1) = datetime('2019-05-27');
S.rat(n+1) = 'K265';
S.tethered(n+1) = true;

%T181
OE = PB_import_opto_ephys_log('T181');
s = struct;
s.date = OE.T181.Date;
s.rat = repmat("T181", numel(s.date),1);
s.tethered = true(numel(s.date),1);
S = structcat(S,s);

% T212
OE = PB_import_opto_ephys_log('T212');
vl_sess = OE.T212.ephys == 1 & OE.T212.poor_performance ~= 1 & OE.T212.few_trials ~= 1;
s = struct;
s.date = OE.T212.Date(vl_sess);
s.rat = repmat('T212', numel(s.date),1);
s.tethered = true(numel(s.date),1);
S = structcat(S,s);

% A230, tethered
s = struct;
A230_in_physrig = bdata('select sessiondate from sessions where sessiondate>="2019-07-15" and ratname="A230" and hostname="Rig205"');
s.date = datetime(A230_in_physrig);
s.rat = repmat('A230', numel(s.date),1);
s.tethered = true(numel(s.date),1);
S = structcat(S,s);

% A242
A242_in_physrig = bdata('select sessid from sessions where sessiondate>="2019-06-18" and ratname="A242" and hostname="Rig205"');
A242_in_physrig_but_maybe_not_tethered = [704440 704117 705282 709379]; % for these sessions there is no physdata collected and video shows rat not plugged in or video is broken so best to exclude these sessions
A242_in_physrig = setdiff(A242_in_physrig, A242_in_physrig_but_maybe_not_tethered);
sessstr='';
for sx=1:numel(A242_in_physrig)
    sessstr=[sessstr, ',' num2str(A242_in_physrig(sx))];
end
s = struct;
s.date = bdata(['select sessiondate from sessions where sessid in (' sessstr(2:end) ')']);
s.date = datetime(s.date);
s.rat = repmat('A242', numel(s.date),1);
s.tethered = true(numel(s.date),1);
S = structcat(S,s);
%% Get the pre-surgery sessions
Surgery_dates = struct;
Surgery_dates.A230 = '2019-07-01';
Surgery_dates.A242 = '2019-06-01';
Surgery_dates.K265 = '2018-05-16';
Surgery_dates.T173 = '2018-01-26';
Surgery_dates.T176 = '2018-04-01';
Surgery_dates.T181 = '2018-05-28';
Surgery_dates.T212 = '2018-08-04';

the_rats = unique(fieldnames(Surgery_dates));
for r = the_rats(:)'; r = r{:};
    date_30d_b4 = datestr(datetime(Surgery_dates.(r))-30, 'yyyy-mm-dd');
    s = struct;
    s.date = bdata(['select sessiondate from sessions where sessiondate >="' date_30d_b4 ...
                    '" and sessiondate <= "' Surgery_dates.(r) '" and ratname="' r '"']);
    s.rat = repmat(r, numel(s.date),1);
    s.tethered = false(numel(s.date),1);
    S = structcat(S, s);
end
%% Get data for comparing n_trials_done and prct_correct
the_rats = unique(S.rat);
Psych = struct;
trials_done = {};
trials_hit = {};
sens = [];
bias = [];
lapse = [];
k = 0;
for r = the_rats(:)'; r = r{:};
    k = k + 1;
    for tethered = [0,1]
        [Psych.(r), ~, Trials] = PB_plot_performance(r, S.date(S.rat == r & S.tethered == tethered), ...
                                                          'suppress_plots', true, ...
                                                          'control_only', true);
        sessids = unique(Trials.sessid);
        for i = 1:numel(sessids)
            vl_done = Trials.sessid == sessids(i) & ...
                      ~Trials.violated & ...
                      Trials.responded & ...
                      Trials.trial_type == 'a';
            vl_hit = vl_done & Trials.is_hit;
            trials_done{k,tethered+1}(i,1) = sum(vl_done);
            trials_hit{k, tethered+1}(i,1) = sum(vl_hit) / sum(vl_done) * 100;
        end
        sens(k, tethered+1) = Psych.(r).sens.est(1);
        bias(k, tethered+1) = Psych.(r).bias.est(1);
        lapse(k, tethered+1) = (Psych.(r).gamma0.est(1)+(100-Psych.(r).gamma1.est(1)))/2;
    end
end
%% Plot parameters
plots_fldr = ['C:/ratter/Analysis/tzluo/Plots/' script_name];
if ~exist(plots_fldr, 'dir'), mkdir(plots_fldr), end
show_legend = false;
show_pval = false;
fig_height = 330; % enforce the same height for all figures
%% Psychometric curves
figure
fig_prepare_axes
set(gcf, 'Position', [800 800, 360, fig_height])
k = 0;
x = -40:40;
for r = the_rats(:)'; r = r{:};
    k = k + 1;
    y = logistic4(x, Psych.(r).prct.gamma0(1), Psych.(r).prct.gamma1(1), Psych.(r).prct.sens(1), Psych.(r).prct.bias(1));
    handles(k) = plot(x,y, 'LineWidth', 1.5);
end
plot(xlim, max(ylim)*[1,1], 'k--', 'LineWidth', 0.5)
plot(xlim, 50*[1,1], 'k--', 'LineWidth', 0.5)
plot(0*[1,1], ylim, 'k--', 'LineWidth', 0.5)
xlabel('#R - L clicks')
ylabel('Fraction chose right')
set(gca, 'XLim', [min(x), max(x)])
if show_legend
    legend(handles, the_rats, 'Location', 'Best')
end
saveas(gcf, [plots_fldr filesep 'behav_performance.png'])
saveas(gcf, [plots_fldr filesep 'behav_performance.svg'])
%% Compare n_trials_done and prct_correct
trials_done_avg = cellfun(@mean, trials_done);
trials_done_sem = cellfun(@sem, trials_done);

figure
fig_prepare_axes
set(gcf, 'Position', [800 800, 360, fig_height])
set(gca, 'DataAspectRatio', [1,1,1])
for r = 1:numel(the_rats)
    hdls(r) = errorbar(trials_done_avg(r,1), trials_done_avg(r,2), ...
             trials_done_sem(r,2),  trials_done_sem(r,2), ...
              trials_done_sem(r,1),  trials_done_sem(r,1), 'capsize', 2, 'linewidth', 1);
end
fig_yMin0
fig_xMin0
quickplot('diag')
xlabel('Untethered')
ylabel('Tethered')
title('Trials completed')
if show_legend, legend(hdls, the_rats, 'Location', 'Best'), end
saveas(gcf, [plots_fldr filesep 'n_trials_completed.png'])
saveas(gcf, [plots_fldr filesep 'n_trials_completed.svg'])
%% Percent correct
prct_hit_avg = cellfun(@nanmean, trials_hit);
prct_hit_sem = cellfun(@sem, trials_hit);

figure
fig_prepare_axes
set(gcf, 'Position', [800 800, 360, fig_height])
set(gca, 'DataAspectRatio', [1,1,1])
for r = 1:numel(the_rats)
    hdls(r) = errorbar(prct_hit_avg(r,1), prct_hit_avg(r,2), ...
                        prct_hit_sem(r,2),  prct_hit_sem(r,2), ...
                        prct_hit_sem(r,1),  prct_hit_sem(r,1), ...
                       '-', 'linewidth', 1, 'capsize', 2);
end
quickplot('diag')
xlabel('Untethered')
ylabel('Tethered')
title('Percent correct')
if show_legend, legend(hdls, the_rats, 'Location', 'Best'), end
saveas(gcf, [plots_fldr filesep 'prct_correct.png'])
saveas(gcf, [plots_fldr filesep 'prct_correct.svg'])
%% Psychometric slope
pval = pval2str(signrank(sens(:,1), sens(:,2)));
figure
fig_prepare_axes
set(gcf, 'Position', [800 800, 360, fig_height])
set(gca, 'DataAspectRatio', [1,1,1])
hdls = [];
for r = 1:numel(the_rats)
    hdls(r) = errorbar(sens(r,1), sens(r,2), ...
                       1.96*sens_std(r,2), 1.96*sens_std(r,2), ...
                       1.96*sens_std(r,1), 1.96*sens_std(r,1), ...
                       '-', 'linewidth', 1, 'capsize', 2);
end
fig_yMin0
fig_xMin0
quickplot('diag')
xlabel('Untethered')
ylabel('Tethered')
if show_pval
    title(['slope (p = ' pval ')'])
else
    title('slope')
end
if show_legend, legend(hdls, the_rats, 'Location', 'Best'), end
saveas(gcf, [plots_fldr filesep 'slope.svg'])
%% Psychometric bias
pval = pval2str(signrank(abs(bias(:,1)), abs(bias(:,2))));
figure
fig_prepare_axes
set(gcf, 'Position', [800 800, 360, fig_height])
set(gca, 'DataAspectRatio', [1,1,1])
hdls = [];
for r = 1:numel(the_rats)
    hdls(r) = plot(abs(bias(r,1))*[1,1], abs(bias(r,2))*[1,1], 'o');
    set(hdls(r), 'MarkerFaceColor', get(hdls(r), 'Color'))
end
% fig_yMin0
% fig_xMin0
quickplot('diag')
xlabel('Untethered')
ylabel('Tethered')
if show_pval
    title(['bias (p = ' pval ')'])
else
    title('bias')
end

if show_legend, legend(hdls, the_rats, 'Location', 'Best'), end
saveas(gcf, [plots_fldr filesep 'bias.svg'])
%% Psychometric lapse
pval = pval2str(signrank(lapse(:,1), lapse(:,2)));
figure
fig_prepare_axes
set(gcf, 'Position', [800 800, 360, fig_height])
set(gca, 'DataAspectRatio', [1,1,1])
hdls = [];
for r = 1:numel(the_rats)
    hdls(r) = errorbar(lapse(r,1), lapse(r,2), ...
                       1.96*lapse_std(r,2), 1.96*lapse_std(r,2), ...
                       1.96*lapse_std(r,1), 1.96*lapse_std(r,1), ...
                       '-', 'linewidth', 1, 'capsize', 2);
end
fig_yMin0
fig_xMin0
quickplot('diag')
xlabel('Untethered')
ylabel('Tethered')
if show_pval
    title(['lapse (p = ' pval ')'])
else
    title('lapse (%)')
end

if show_legend, legend(hdls, the_rats, 'Location', 'Best'), end
saveas(gcf, [plots_fldr filesep 'lapse.svg'])