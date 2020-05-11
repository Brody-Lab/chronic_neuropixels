if ~exist('S', 'var')
    load([P.data_folder_path filesep 'exp_decay_mdl_selection_2020_04_28.mat']);
end
big_figure
n_var = size(S.boot.b,1);
k = 0;
idx = S.b(7,:) < 0.2;

idx_reg = 1:(numel(P.ED_trode_regressors_all));
for i = 1:n_var
for j = 1:n_var
    k = k + 1;
    subplot(n_var, n_var, k);
    hold on
    plot(S.boot.b(i,:),S.boot.b(j,:),'ko')
    fig_yMin0
    fig_xMin0
    if i == 1
        title(P.ED_trode_regressors_all{k})
    end
end
end
%%
big_figure
k = 0;

for i = 1:n_var
    k = k + 1;
    subplot(5,5,k)
    set(gca, 'nextplot', 'add')
    for j = 1:size(S.boot.cil,2)
        plot(j*[1,1], [S.boot.cil(i,j), S.boot.ciu(i,j)], 'k-')
        plot(j, S.boot.b(i,j), 'ro')
    end
    title(varnames{k})
    refline(0,0)
    
    shadeplot(xlim, S.cil(i)*[1,1], S.ciu(i)*[1,1], 'Color', [0,0,1])
    shadeplot(xlim, cil(i)*[1,1], ciu(i)*[1,1], 'Color', [1,0,0])
    plot(xlim, S.b(i)*[1,1], 'm')
end