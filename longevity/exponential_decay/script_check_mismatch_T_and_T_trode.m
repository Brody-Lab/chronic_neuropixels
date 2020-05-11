P=get_parameters;
T_trode = make_T_trode(Cells);
T = get_metrics_from_Cells(Cells);
%%
sum(T_trode.unit)/size(T_trode,1)
sum(T.unit)/sum(T.n_elec)
%%
for i = 1:size(T,1)
    idx = T_trode.identifier == T.identifier(i) & ...
          T_trode.days_elapsed == T.days_elapsed(i);
    if sum(idx) ~= T.n_elec(i)
        error('electrode # mismatch %i', i)
    end
    if sum(T_trode.unit(idx)) ~= T.unit(i)
        error('unit # mismatch %i', i)
    end
end

unique_id = unique(T_trode(:, {'identifier', 'days_elapsed'}), 'rows');
for i = 1:size(unique_id,1)
    idx_trode = T_trode.identifier == unique_id.identifier(i) & ...
                T_trode.days_elapsed == unique_id.days_elapsed(i);
    idx_sess = T.identifier == unique_id.identifier(i) & ...
                T.days_elapsed == unique_id.days_elapsed(i);
    if sum(idx_trode) ~= T.n_elec(idx_sess)
        error('electrode # mismatch %i', i)
    end
    if sum(T_trode.unit(idx_trode)) ~= T.unit(idx_sess)
        error('unit # mismatch %i', i)
    end
end
%%
plot_ED_trode_fit(S);
plot_average_stability(T, 'metric', 'unit', ...
                              'axes', gca, ...
                            'normalize_by_electrodes', true);
ylim([0,2])