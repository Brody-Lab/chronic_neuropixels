



function model_implantation()

analyze_gain_noise_data;

bank_0_noise_new = bank_0_noise{idx_new};
bank_1_noise_new = bank_1_noise{idx_new};
bank_0_noise_exp = bank_0_noise{~idx_new};
bank_1_noise_exp = bank_1_noise{~idx_new};


min_val=1e-1;
noise(1).exp=bank_0_noise_exp(bank_0_noise_exp>min_val);
noise(2).exp=bank_1_noise_exp(bank_1_noise_exp>min_val);
noise(1).new=bank_0_noise_new(bank_0_noise_new>min_val);
noise(2).new=bank_1_noise_new(bank_1_noise_new>min_val);

mu = [1 1]*mean(log10([noise(1).new(:);noise(2).new(:)]));
sigma = cov(log10(noise(1).new(:)),log10(noise(2).new(:)));
sigma(1,1) = var(log10([noise(1).new(:);noise(2).new(:)]));
sigma(2,2) = var(log10([noise(1).new(:);noise(2).new(:)]));
sigma(1,2)=0.0074;
sigma(2,1)=0.0074;

x=logspace(-3,-0,20);
count=0;
n_samples=1e5;
clear M r2
for i=x(:)'
    count=count+1;    
    for k=1:50
        predicted_noise = mvnrnd(mu,sigma,n_samples) + i.*mvnrnd(mu*i,[sigma(1) 0; 0 sigma(4)],n_samples);
        r2(count,k) = rsquare(zscore(10.^predicted_noise(:,1)),zscore(10.^(predicted_noise(:,2))));
        M(count,k) = median(10.^predicted_noise(:));
    end
end
a=1;
end
