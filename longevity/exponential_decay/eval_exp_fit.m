function [yhat] = eval_exp_fit(b, T_trode, t)
    b=b(:);
    [X0, Xt] = make_X0_Xt(T_trode, ~isnan(b));
    b = b(~isnan(b));
    yhat = eval_exp(b,X0,Xt,t);
end

function yhat = eval_exp(b,X0,Xt,t)
    n = size(X0,2);
    y0 = X0*b(1:n);
    tDtau = t ./ (Xt*b(n+1:end));
    yhat = y0.*exp(-tDtau);   
end