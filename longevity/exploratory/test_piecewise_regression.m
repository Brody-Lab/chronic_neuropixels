T = get_metrics_from_Cells(Cells, 'condition_on', {'AP'});
T = T(T.days_elapsed > 0,:);
T.metric = T.unit./T.n_elec;
%%
figure
T1 = T(T.condition==1,:);
plot(T1.days_elapsed+1, T1.metric, 'o')
hold on
T2 = T(T.condition==2 & T.days_elapsed<300,:);
plot(T2.days_elapsed+1, T2.metric, 'o')
set(gca, 'xscale', 'log')
%%
x = log(T1.days_elapsed);
y = log(T1.metric);

slm = slmengine(x, y, ...
                'knots', [0 6], ...
                'degree',1,...
                'plot','off', ...
                'interiorknots', 'fixed');
pp = slm2pp(slm);
idx = T1.days_elapsed >= slm.knots(2);
x2 = x(idx);
ssx = sum(x2.^2) - sum(x2)^2/numel(x2);
yhat = slmeval(x2, slm);
y2 = y(idx);
mse = sum((y2-yhat).^2)/(sum(idx)-2);
se = sqrt(mse/ssx);
tcv = tinv(0.975, sum(idx)-2); % T-CDF's critical value
err95 = se*tcv;
fprintf('\n[%f, %f]', -1/(pp.coefs(2)-err95), -1/(pp.coefs(2)+err95))

sum((y-slmeval(x,slm)).^2)

figure
plot(exp(x), exp(y), 'ko')
hold on
plot(exp(sort(x)), exp(slmeval(sort(x), slm)))
set(gca, 'xscale', 'log','yscale', 'log')
%%
x = log(T2.days_elapsed);
y = log(T2.metric);

slm = slmengine(x, y, ...
                'knots', [0 6], ...
                'degree',1,...
                'plot','off', ...
                'interiorknots', 'fixed');
pp = slm2pp(slm);
idx = T1.days_elapsed >= slm.knots(2);
x2 = x(idx);
ssx = sum(x2.^2) - sum(x2)^2/numel(x2);
yhat = slmeval(x2, slm);
y2 = y(idx);
mse = sum((y2-yhat).^2)/(sum(idx)-2);
se = sqrt(mse/ssx);
tcv = tinv(0.975, sum(idx)-2); % T-CDF's critical value
err95 = se*tcv;
fprintf('\n[%f, %f]', -1/(pp.coefs(2)-err95), -1/(pp.coefs(2)+err95))

sum((y-slmeval(x,slm)).^2)

figure
plot(exp(x), exp(y), 'ko')
hold on
plot(exp(sort(x)), exp(slmeval(sort(x), slm)))
% set(gca, 'xscale', 'log','yscale', 'log')

%%
x = log(T2.days_elapsed);
y = log(T2.metric);
[x, idx] = sort(x);
y=y(idx);

m =[x,ones(numel(x),1)]\y;
yhat = [x,ones(numel(x),1)]*m;
n = sum(idx);
mse = sum((y-yhat).^2)/(n-2);
se = sqrt(mse/ssx);
tcv = tinv(0.975, n-2); % T-CDF's critical value
err95 = se*tcv;
% fprintf('\tau: [%f, %f]\n', -1/(pp.coefs(2)-err95), -1/(pp.coefs(2)+err95))
fprintf('\nExponential fit, MSE=%f\n', mse)

figure
plot(exp(x), exp(y), 'ko')
hold on
plot(exp(x), exp(yhat))
set(gca, 'xscale', 'log','yscale', 'log')
%% Exponential fit
x = T2.days_elapsed;
y = log(T2.metric);
[x, idx] = sort(x);
y=y(idx);

m =[x,ones(numel(x),1)]\y;
yhat = [x,ones(numel(x),1)]*m;
n = sum(idx);
mse = sum((y-yhat).^2)/(n-2);
se = sqrt(mse/ssx);
tcv = tinv(0.975, n-2); % T-CDF's critical value
err95 = se*tcv;
% fprintf('\tau: [%f, %f]\n', -1/(pp.coefs(2)-err95), -1/(pp.coefs(2)+err95))
fprintf('\nExponential fit, MSE=%f\n', mse)
figure
plot(x, exp(y), 'ko')
hold on
plot(x, exp(yhat))
set(gca, 'yscale', 'log')