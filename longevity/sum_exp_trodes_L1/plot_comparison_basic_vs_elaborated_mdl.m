add_folders_to_path
P=get_parameters;
load(P.sum_exp_trodes.orig_vs_elab_data_path, 'diff_boot', 'diff_obsv')
figure
set(gca, 'NextPlot', 'Add', ...
         'fontsize', 14)
hdl_boot = histogram(diff_boot);
hdl_obsv = plot(diff_obsv*[1,1], ylim);
xlabel('LL(elaborated) - LL(basic)')
ylabel('# bootstrap samples')
legend([hdl_boot, hdl_obsv], {'Boostrapped samples', 'Observed'}, 'location', 'best')
pval = sum(abs(diff_boot)>=abs(diff_obsv))/numel(diff_boot);
if pval == 0
    pval = ['p < ' num2str(1/numel(diff_boot))];
else
    pval = ['p \leq ' num2str(1/numel(diff_boot))];
end
title(pval)