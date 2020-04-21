% EXAMPLE_PROBE_BANK_NOISE_SCATTER make a scatter plot of noise on bank 0
% versus 1 for an example probe
%
%=OPTIONAL
%       axes
%           An axes object. If this were not provided, a new figure is made
function[] = example_probe_bank_noise_scatter(probe_sn,varargin)
parseobj = inputParser;
addParameter(parseobj, 'axes', [], @(x) isempty(x) || isa(x,  'matlab.graphics.axis.Axes'));
parse(parseobj, varargin{:});
P_in = parseobj.Results;
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
analyze_gain_noise_data;
%% plot
if isempty(P_in.axes)
    figure('Position', P.figure_position_gn_summary)
else
    axes(P_in.axes)
end
idx=find(cell2mat(cellfun(@(x)x==P_in.probe_sn,probe_sn,'uni',0)));
A = bank_0_noise{idx};
B = bank_1_noise{idx};
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