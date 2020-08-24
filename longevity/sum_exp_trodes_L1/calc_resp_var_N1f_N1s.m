function yhat = calc_resp_var_N1f_N1s(b, XN1f, XN1s, Xk, t)
% CALC_RESP_VAR calculate the response variable using the model that has
% separate terms for N1f and N1s
%
%=INPUT
%
%   b
%       A vector of parameter coefficients
%
%   XN1f
%       Design matrix for the N1f term
%
%   XN1s
%       Design matrix for the N1s term
%
%   Xk
%       Design matrix for the k term
%
%   t
%       Number of days since implant minus one (t=num_days-1)
%
%=OUTPUT
%
%   yhat
%       predicted response variable

    n_XN1f = size(XN1f,2);
    n_XN1s = size(XN1s,2);
    kf = b(1);
    ks = b(2);
    bN1f = b(3:n_XN1f+2);
    bN1s = b(n_XN1f+3:n_XN1f+n_XN1s+2);
    bk   = b(n_XN1f+n_XN1s+3:end);
    
    yhat = (XN1f*bN1f.*exp(kf*t) + XN1s*bN1s.*exp(ks*t)) .* exp(t.*Xk*bk); 
end