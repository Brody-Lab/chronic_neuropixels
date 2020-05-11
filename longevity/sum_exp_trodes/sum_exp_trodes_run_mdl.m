% SUM_EXP_TRODES_RUN_MDL A scrip for running the sum-of-exponentials model
% on the unit count data from electrodes

tic
P = get_parameters;
if ~exist('Cells', 'var')
    load(P.Cells_path)
end
S = sum_exp_trodes_select_mdl(Cells, 'Iterations', 10, 'KFold', 5);
S = sum_exp_trodes_compute_CI(S, 'i_mdl', 1:10);
save(P.sum_exp_trodes.data_path, 'S')
fprintf('\nDone after %0.f seconds!\n', toc)