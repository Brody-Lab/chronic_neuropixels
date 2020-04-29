% BANK_NOISE_SCATTER make a scatter plot of noise on bank 0
% versus 1 for all probes or just one example probe
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = bank_noise_scatter(varargin)
parseobj = inputParser;
addParameter(parseobj,'probe_sn',[]);
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
if ~isempty(P_in.probe_sn)
    P_in.probe_sn = probe_sn; clear probe_sn;
    P=get_parameters;
    T= readtable(P.gain_noise_log_path);
    idx = T.probe_sn == P_in.probe_sn;
    ind = find(idx);
    if isempty(ind)
        error('Probe S/N not in gain noise log.');
    end
    data_file_path = [P.gain_noise_fldr_path filesep T.recording_id{ind} '.csv'];
    D = readtable(data_file_path);
    idx_channels = ((1:384) + 384)<=T.electrodes_implanted(ind);
    min_channels=50;
    if sum(idx_channels)<min_channels                 
        error('Bank 1 wasn''t in the brain for this probe.');
    end
end
analyze_gain_noise_data;
%% plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
if ~isempty(P_in.probe_sn)
    idx=find(cell2mat(cellfun(@(x)x==P_in.probe_sn,probe_sn,'uni',0)));
    A = bank_0_noise{idx};
    B = bank_1_noise{idx};
    r2=rsquare(A,B);    
    scatter(A,B,'.');
    set(gca,P.axes_properties{:},'xscale','log','yscale','log');
    l=refline(1,0);
    h=scatter(A,B,'o');
    h.SizeData=20;
    h.MarkerFaceColor='r';
    h.MarkerEdgeColor=[1 1 1];
    h.LineWidth=0.1;
    l.LineStyle =':';
    xlabel('Bank 0 Noise (\muV_R_M_S)');
    ylabel('Bank 1 Noise (\muV_R_M_S)');    
    text(sprintf('R^2 = %0.3f',rsquare(A,B)));
else
    lims=[-1 2];
    A=cat(1,bank_0_noise_z{:});
    B=cat(1,bank_1_noise_z{:});
    r2=rsquare(A,B);    
    A(A<lims(1))=lims(1);
    A(A>lims(2))=lims(2);
    B(B<lims(1))=lims(1);
    B(B>lims(2))=lims(2);    
    scatter(A,B,7,'ok','MarkerFaceAlpha',0.2,'MarkerFaceColor','k','MarkerEdgeColor','none');hold on;
    set(gca,P.axes_properties{:},'xscale','linear',...
        'yscale','linear','ylim',lims,'xlim',lims,...
        'ytick',lims(1):lims(2),'xtick',lims(1):lims(2),...
        'yticklabel',{'\leq-1', '0','1' '\geq2'},...
        'xticklabel',{'\leq-1', '0','1' '\geq2'});
    l=refline(1,0);    l.LineStyle =':';l.LineWidth=1;l.Color='k';
    xlabel({'Bank 0 RMS Noise','z-scored (AU)'});
    ylabel({'Bank 1 RMS Noise','z-scored (AU)'}); 
    text(-0.8,1.7,sprintf('R^2 = %0.2f',r2));    
end


