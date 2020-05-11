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