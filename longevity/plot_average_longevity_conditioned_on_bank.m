% plot_average_longevity_conditioned_on_bank
% plot total number of recorded SUs over time across sessions, with temporal binning, separately for banks 1 and 0.
P = get_parameters;
days_elapsed = cellfun(@(x)x.days_since_surgery,Cells,'UniformOutput',1)+1;
n_good_units = cellfun(@(x)x.n_good_units,Cells,'UniformOutput',1);    
unique_banks=cellfun(@(x) x.unique_bank, Cells);
bin_edges = P.longevity_time_bin_edges;
bin_centers = bin_edges(1:end-1)+diff(bin_edges)/2;
n_time_bins = numel(bin_edges)-1;
figure;
for k=0:1
    clear boots
    in_bank = unique_banks==k;
    for i=1:n_time_bins
        idx = days_elapsed>bin_edges(i) & days_elapsed<bin_edges(i+1);
        boots(:,i) = bootstrp(1000,@mean,n_good_units(idx(:) &  in_bank(:)   ));
    end
    t(k+1)=shadedErrorBar(bin_centers,boots,{@mean,@std});hold on;
    t(k+1).mainLine.LineWidth=2;
    t(k+1).mainLine.Color = bank_colors(k+1,:);
    t(k+1).patch.FaceAlpha=0.5;
    t(k+1).patch.FaceColor = bank_colors(k+1,:);
    set(gca,'xscale','log',...
            'yscale','log',...
            'xlim',[1,mean(bin_edges(end-1:end))],'xtick',[1 10 50 100 200 300 400],'ylim',[50 400]);
    xlabel('Days Since Implant');ylabel('No. SU');        
end
legend([t.mainLine],{'Bank 0','Bank 1'});