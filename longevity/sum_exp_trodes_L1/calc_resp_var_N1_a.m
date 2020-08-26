function yhat = calc_resp_var_N1_a(b, XN1, Xk, t)
% CALC_RESP_VAR calculate the response variable using the model that has
% separate terms for N1 and a
%
%=INPUT
%
%   b
%       A vector of parameter coefficients
%
%   XN1
%       Design matrix for the N1 term
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

    n_XN1 = size(XN1,2);
    kf = b(1);
    ks = b(2);
    a = (cos(b(3))+1)/2;
    bN1 = b(4:n_XN1+3);
    bk = b(n_XN1+4:end);
    
    N1 = XN1*bN1;
    
    yhat = (a*N1.*exp(kf*t) + (1-a)*N1.*exp(ks*t)) .* exp(t.*Xk*bk); 
end