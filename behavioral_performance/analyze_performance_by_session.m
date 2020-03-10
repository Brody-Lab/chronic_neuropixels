% ANALYZE_PERFORMANCE_BY_SESSION provides summary statistics for each
% session in BEHAVIOR_TABLE.CSV
function [] = analyze_performance_by_session()
    add_folders_to_path;
    P = get_parameters;
    T=readtable(P.behavior_table_path);    
    [~,idx] = sort(T.sessid);
    T = T(idx,:);
    sessids = concatenate_for_sql(T.sessid);
    fprintf('Fetching protocol data...')
    PD=bdata(['select protocol_data from sessions s where s.sessid in (' sessids ')']);
    fprintf(' done\n')
    for i = 1:numel(PD)
        PD{i} = remove_trials_from_pd(PD{i});
    end
    
    B = struct;
    B.trials_done = cellfun(@(x) numel(x.hits), PD);
    B.prct_hit = cellfun(@(x) sum(x.hits)/numel(x.hits), PD)*100;
    B.sens = nan(numel(PD,1));
    B.bias = nan(numel(PD,1));
    B.lapse = nan(numel(PD,1));
    for i = 1:numel(PD)
        PD{i}.pokedR = (PD{i}.sides == 'r' &  PD{i}.hits) | ...
                       (PD{i}.sides == 'l' & ~PD{i}.hits);
        Psych = fit_psychometric(PD{i});
        B.sens(i,1) = Psych.sens;
        B.bias(i,1) = abs(Psych.bias);
        B.lapse(i,1) = Psych.lapse;
    end
    B = struct2table(B);
    B = [T, B];
    writetable(B, P.performance_by_session_path);
end
%% remove_trials
% Ensure all fields of the protocol data structure to have the same number
% of trials and remove trials with violations or optogenetic illumination
%=INPUT
%
%   pd
%       a structure with protocol data for one session
%
%=OUTPUT
%
%   pd
%       a structure with protocol data for one session
function pd = remove_trials(pd)
    n_trials = min(structfun(@numel, pd));
    pd = structfun(@(x) x(1:n_trials), pd, 'uni', 0);
    is_on = cellfun(@(x) x.ison, pd.stimdata);
    pd = structfun(@(x) x(~is_on), pd, 'uni', 0);
    pd = structfun(@(x) x(~pd.violations), pd, 'uni', 0);
    clicks_gamma=cellfun(@(x) x.gamma, pd.bupsdata);
    requires_accumulation = abs(clicks_gamma) < 5;
    pd = structfun(@(x) x(requires_accumulation), pd, 'uni', 0);
end