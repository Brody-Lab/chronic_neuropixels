function [] = debug_exp_fit(b,X0train,Xttrain,ttrain, ytrain)
close all
t_unique = unique(ttrain); 
y_train_hat= eval_exp(b,X0train,Xttrain,ttrain);
for i =1:numel(t_unique)
    data_avg(i)=mean(ytrain(ttrain==t_unique(i))); 
    hat_avg(i)=mean(y_train_hat(ttrain==t_unique(i))); 
end 
figure
set(gca, 'xscale', 'log')
hold on
plot(t_unique, data_avg,'ko'); 
plot(t_unique, hat_avg,'ro'); 
end
%% Evaluate a exponential function
function yhat = eval_exp(b,X0,Xt,t)
    n = size(X0,2);
    y0 = X0*b(1:n);
    tDtau = t ./ (Xt*b(n+1:end));
    yhat = y0.*exp(-tDtau);   
end