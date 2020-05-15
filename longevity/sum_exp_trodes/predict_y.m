% CALC_LAMBDA calculate lambda for the sum-of-exponentials model
%
%=INPUT
%
%   XN1
%       Design matrix in the N1 term
%
%   Xk
%
%       Design matrix in the decay rate term
%
%
%   t
%
%       days from initial day
%
%   y
%
%       the response variable (unit count, etc.)
function yhat = predict_y(b, XN1, Xk, t)
    n_XN1 = size(XN1,2);
    a = b(1);
    kf = b(2);
    ks = b(3);
    bN1 = b(4:n_XN1+3);
    bk = b(n_XN1+4:end);
    N1 = XN1*bN1;
    yhat = N1 .* (a*exp(kf*t) + (1-a)*exp(ks*t)) .* (exp(t.*Xk*bk) ); 
end