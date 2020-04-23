T = get_metrics_from_Cells(Cells, 'condition_on', {'DV', 'EI', 'AP'});
metric = 'unit';
T.metric= T.(metric) ./ T.n_elec;

T2 = T(:, {'days_elapsed', 'trode_DV_avg', 'trode_AP_avg', 'trode_EI_avg', metric});

T2.days_elapsed=log(T2.days_elapsed)

eqn = [metric ' ~ days_elapsed + days_elapsed*trode_DV_avg + ' ...
                 'days_elapsed*trode_AP_avg + days_elapsed*trode_EI_avg'];

fitlm(T2, eqn)

%%
T = get_metrics_from_Cells(Cells, 'condition_on', {'DV', 'EI'});
metric = 'unit';
T.metric= T.(metric) ./ T.n_elec;


T2 = T(:, {'days_elapsed', 'trode_DV_avg', 'trode_EI_avg', metric});

eqn = [metric ' ~ days_elapsed + days_elapsed:trode_DV_avg + ' ...
                 ' days_elapsed:trode_EI_avg'];
fitlm(T2, eqn)

%%
T_trode = struct;
k = 0;
for i =1:numel(Cells)
    n_trode = sum(Cells{i}.electrodes.in_brain);
    T_trode.days_elapsed{i,1} = repmat(Cells{i}.days_since_surgery, n_trode, 1);
    T_trode.DV{i,1} = Cells{i}.electrodes.DV;
    T_trode.AP{i,1} = Cells{i}.electrodes.AP;
    T_trode.ML{i,1} = Cells{i}.electrodes.ML;
    T_trode.shank_pos{i,1} = ceil(Cells{i}.electrodes.index/2);
    for j = 1:n_trode
        T_trode.unit{i,1}(j,1) = sum(Cells{i}.electrode == Cells{i}.electrodes.index(j));
        T_trode.single_unit{i,1}(j,1) = sum(Cells{i}.electrode == Cells{i}.electrodes.index(j) & ...
                                     Cells{i}.ks_good(:));
    end
end
T_trode = structfun(@cell2mat, T_trode, 'uni', 0);
T_trode = struct2table(T_trode);
%%
metric = 'single_unit';
eqn = [metric ' ~ days_elapsed*DV ' ...
               '+ days_elapsed*AP' ...
               '+ days_elapsed*shank_pos' ...
                ];
fitlm(T_trode, eqn)