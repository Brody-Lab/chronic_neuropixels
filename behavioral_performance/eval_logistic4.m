function [y] = eval_logistic4(x, gamma0, gamma1, sens, bias)
    y = gamma0 + gamma1./(1 + exp(-sens*(x-bias)));