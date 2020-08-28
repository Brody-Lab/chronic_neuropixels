add_folders_to_path % add the "chronic_neuropixels" to MATLAB search path

% panel A
show_equations

% panel B
plot_elab_mdl_coeff

% text
[prct_larger, pval] = report_on_basic_vs_elaborated_mdl;
str = sprintf('The log-likelihood of the elaborated model is higher than the basic model %0.1f%% (p < %0.3f)', ...
              prct_larger, pval);
msgbox(str)