% FIT_PSYCHOMETRIC fit a 4-parameter logistic functon
%
%=INPUT
%
%   pd
%       a structure with protocol data for one session
%
%=OUTPUT
%
%   Psych
%       a structure with information on the fitted psychometric
function Psych = fit_psychometric(pd)
    Psych = fit_logistic4(pd.pokedR,pd.n_right-pd.n_left);
    Psych.gamma0=Psych.gamma0*100;
    Psych.gamma1=Psych.gamma1*100;
    lapse1 = Psych.gamma0;
    lapse2 = 100-(Psych.gamma0+Psych.gamma1);
    Psych.lapse = (lapse1+lapse2)/2;
end