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
function pd = remove_trials_from_pd(pd)
    n_trials = min(structfun(@numel, pd));
    pd = structfun(@(x) x(1:n_trials), pd, 'uni', 0);
    is_on = cellfun(@(x) x.ison, pd.stimdata);
    pd = structfun(@(x) x(~is_on), pd, 'uni', 0);
    pd = structfun(@(x) x(~pd.violations), pd, 'uni', 0);
    clicks_gamma=cellfun(@(x) x.gamma, pd.bupsdata);
    requires_accumulation = abs(clicks_gamma) < 5;
    pd = structfun(@(x) x(requires_accumulation), pd, 'uni', 0);
end