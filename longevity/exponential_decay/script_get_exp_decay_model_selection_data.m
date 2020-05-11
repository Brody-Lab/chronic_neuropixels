tic
S = select_exp_mdl(Cells, 'Iterations', 25, 'KFold', 10);
S = est_coeff_of_sel_mdl(S, 'n_boot', 1000);
toc;
save([P.data_folder_path filesep 'exp_decay_mdl_selection.mat'], 'S')

