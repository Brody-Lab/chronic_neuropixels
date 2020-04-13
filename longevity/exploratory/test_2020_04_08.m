% generate toy data
x = 2.^(0:0.1:10);
p = [5,4,0,2000];
sigma2 = 1;
n = numel(x);
y = sum_2_exp_decay(x,p) + randn([1,n])*sqrt(sigma2);
y(y<0)=0;
figure
plot(x,y, 'o')
set(gca, 'xscale', 'log')
%%
p0 = [1,1,1,100];
negLL_sum_2_exp_decay(x,y,p0);

lower_bounds = [0,0,0,0];
upper_bounds = inf*ones(1,4);

[p_hat, p_CI, LL, hessian] = fit_sum_2_exp_decay(x,y, 'upper_bounds', [inf,inf,0,0]);

hold on
plot(x, sum_2_exp_decay(x,p_hat))
for i = 1:4
    fprintf('\n generative p(%i) = %0.2f; estimated p(%i) = %0.2f [%0.2f, %0.2f]', ...
            i, p(i), i, p_hat(i), p_CI(i,1), p_CI(i,2))
end