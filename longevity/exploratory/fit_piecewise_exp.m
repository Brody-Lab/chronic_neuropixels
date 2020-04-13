function S = fit_piecewise_exp(x,y)

x=log(x);
y=log(y);
x_break0 = log(4);

[x_fmin, f_fmin0] = fminbnd(@(x_break) fit_piecewise_lin(x, y, x_break), min(x), max(x))
[~, p0] = fit_piecewise_lin(x, y, x_fmin);

[x_fmin, f_fmin1, ~, ~, ~, ~, hessian] = ...
        fmincon(@(x_break) fit_piecewise_lin(x, y, x_break), ...
                x_break0, ...
                [], [], [], [], ...
                min(x), log(7), [])




lb = [min(x),  0, -inf, -inf];
ub = [max(x),  max(y), inf, inf];

% lb = p0-p0/10;
% ub = p0+p0/10;

[x_fmin, f_fmin2, exitflag, output, ~, grad, hessian] = ...
        fmincon(@(p) fit_piecewise_lin_all(x, y, p), ...
                p0, ...
                [], [], [], [], ...
                lb, ub, [], ...
                optimset('MaxIter', 200, ...
                         'Algorithm', 'interior-point', ...
                         'Display', 'off'));
f_fmin0
f_fmin1
S.x_fmin = exp(x_fmin);
SE = sqrt(diag(hessian^-1));
S.CI = S.x_fmin + SE * 1.96*[-1, 1];
S.CI = exp(S.CI);
S.name = {'break', 'intercept', 'slope1', 'slope2'};
end

function [sse, m] = fit_piecewise_lin(x,y,x_break)
idx=x<x_break;
x1=x(idx);
y1=y(idx);
x2=x(~idx);
y2=y(~idx);

m1 =[x1,ones(sum(idx),1)]\y1;
m2 =[x2,ones(sum(~idx),1)]\y2;

y1_hat = [x1,ones(sum(idx),1)]*m1;
y2_hat = [x2,ones(sum(~idx),1)]*m2;

sse = sum((y1_hat-y1).^2) + sum((y2_hat-y2).^2);

m = [x_break; m1(2); m1(1); m2(1)];
end

% p(1) = x_break
% p(2) = y0
% p(3) = m1
% p(4) = m2
function sse = fit_piecewise_lin_all(x,y,p)
idx=x<p(1);
x1=x(idx);
y1=y(idx);
x2=x(~idx);
y2=y(~idx);

y1_hat = p(2) + p(3)*x1;
y2_hat = p(2) + p(3)*p(1) + p(4)*x2;

sse = sum((y1_hat-y1).^2) + sum((y2_hat-y2).^2);
end
