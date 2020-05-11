% MAKE_X_XT design matrices for the exponential decay fits
%
%=INPUT
%
%   T_trode
%       Table made using MAKE_T_TRODE
%
%   regressors_used
%       A logical row array specifying whether each regressor is used
%
%   trodes_used
%       A logical vector array specifying whether electrodes are used
%=OUTPUT
%
%   X
%       The columns of the design matrix that specifies the initial count
%
%   Xt
%       The columns of the design matrix that specifies the decay rate,
%       multiplied by t
%
%   tDX
%       The columns of the design matrix that specifies the decay rate, divided by t 
function [X, Xt, tDX] = make_X_Xt(T_trode, regressors_used, trodes_used)
P = get_parameters;
% remove the 1st and n/2+1 column, which are constants
n=numel(regressors_used);
idx_y0 = regressors_used(2:n/2);
idx_k = regressors_used(n/2+2:end);
T_trode = T_trode(trodes_used,:);
X_factors = T_trode{:, P.ED_trode_factors};
X = X_factors(:,idx_y0);
Xt = X_factors(:,idx_k);
X  = [ones(sum(trodes_used),1), X];
Xt = [ones(sum(trodes_used),1), Xt];
tDX = Xt ./ T_trode.days_elapsed;
Xt = Xt .* T_trode.days_elapsed;