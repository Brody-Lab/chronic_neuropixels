% PLOT_AVERAGE_LONGEVITY plot the neuronal stability metrics averaged
% across recording sessions-conditions. Errorbars indicate the standard
% error of the mean.
%
%=INPUT
%
%   Cells
%       A cel array of structures made by COLLECT_CELL_FILES and
%       POSTPROCESS_CELLS
%
%=OPTIONAL INPUT
%
%   anatom_bin_edges
%       If CONDITION_ON were 'ML','AP', or 'DV', then this increasing
%       vector specifies the bin edges for binning anatomical coodinates
%
%   condition_on
%       A char array specifying the different conditions used for
%       averaging. The option are {'bank', 'ML', 'AP', 'DV'}.
%
%   FaceAlpha
%       A scalar value between 0 and 1 that specifies the opacity of the
%       shading that represent the standard error of mean. When a large
%       number of conditions are plotted, this value can be set to 0 for
%       better clarity.
%
%   metric
%       The neuronal stability metric to plot: {'unit', 'single_unit',
%       'event_rate', which is the total firing rate, or 'Vpp',
%       peak-to-peak amplitude of the spike waveform
%
%   normalize_by_electrodes
%       A scalar logical specifying whether to normalize neuronal metrics
%       by the number of electrodes that satisfy the same condition. This
%       is should be turned on CONDITION_ON is not empty
function [] = plot_average_longevity(Cells, varargin)
parseobj = inputParser;
addParameter(parseobj, 'anatom_bin_edges', 0:2:10, ...
    @(x) validateattributes(x, {'numeric'}, {'increasing', 'vector'}))
addParameter(parseobj, 'condition_on', '', @(x) ismember(x, {'bank', ...
                                                             'ML', ...
                                                             'AP', ...
                                                             'DV'}))
addParameter(parseobj, 'metric', 'unit', @(x) ismember(x, {'unit', ...
                                                           'single_unit', ...
                                                           'event_rate', ...
                                                           'Vpp'}))
addParameter(parseobj, 'normalize_by_electrodes', false, @(x) islogical(x)&&isscalar(x))
addParameter(parseobj, 'FaceAlpha', 0.2, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', '>=', 0, '<=', 1}))
parse(parseobj, varargin{:});
P_in = parseobj.Results;
P = get_parameters;
t_bin_edges = P.longevity_time_bin_edges;
t_bin_centers = P.longevity_time_bin_centers;
n_boots = P.longevity_n_boots;
%% Get the appropriate metric
switch P_in.metric
    case 'single_unit'
        Metric = cellfun(@(x) x.ks_good,Cells,'uni',0);
        y_label = 'Single units';
        func = @sum;
    case 'unit'
        Metric = cellfun(@(x) ones(x.n_units,1), Cells, 'uni', 0);
        y_label = 'Units';
        func = @sum;
    case 'event_rate'
        Metric = cellfun(@(x)x.fr,Cells,'uni', 0);
        y_label = 'Total spike rate';
        func = @sum;
    case 'Vpp'
        Metric = cellfun(@(x) x.unitVppRaw,Cells,'uni', 0);
        y_label = 'Peak-to-peak amplitude (uV)';
        func = @mean;
        P_in.normalize_by_electrodes = false;
end
Metric = cellfun(@(x) x(:),Metric,'uni', 0);
%% provide an index for a condition for each cell for each recording
switch P_in.condition_on
    case 'bank'
        for i = 1:numel(Metric)
            I_cells{i,1} = Cells{i}.bank+1;
            I_elec{i,1} = Cells{i}.electrodes.bank+1;
        end
        condition_names = ["Bank 0"; "Bank 1"];
    case {'AP', 'ML', 'DV'}        
        for i = 1:numel(Metric)
            I_cells{i,1} = discretize(Cells{i}.(P_in.condition_on), P_in.anatom_bin_edges);
            if any(isnan(I_cells{i}))
                warning('The bin edges cover only a subset of the range of %s positions.',P_in.condition_on)
                I_cells{i} = I_cells{i}(~isnan(I_cells{i}));
                Metric{i} = Metric{i}(~isnan(I_cells{i}));
            end
            I_elec{i,1} = discretize(Cells{i}.electrodes.(P_in.condition_on), P_in.anatom_bin_edges);
        end
        for i = 1:numel(P_in.anatom_bin_edges)-1
            edge1=num2str(P_in.anatom_bin_edges(i));
            edge2=num2str(P_in.anatom_bin_edges(i+1));
            
            condition_names{i,1} = [P_in.condition_on ' [' edge1 ',' edge2 ']'];
        end        
    otherwise
        I_cells = cellfun(@(x) ones(x.n_units,1), Cells, 'uni', 0);
        I_elec=cellfun(@(x) ones(x.n_electrodes_in_brain,1), Cells, 'uni', 0);
        condition_names="all";
end
I_cells = cellfun(@(x) x(:),I_cells,'uni', 0);
% separate cells within each recording session by its condition
Condition_id=[];
for i = 1:numel(Metric)
    Vals{i,1} = splitapply(func, Metric{i}, findgroups(I_cells{i}));
    unique_bin_inds = unique(I_cells{i});
    Condition_id = [Condition_id; condition_names(unique_bin_inds)];
    Days_elapsed{i,1} = ones(numel(unique_bin_inds),1)*Cells{i}.days_since_surgery;
    for j = 1:numel(unique_bin_inds)
        N_electrodes{i,1}(j,1) = sum(I_elec{i}==unique_bin_inds(j));
        if N_electrodes{i,1}(j,1)==0
            error('Must have an electrode with this condition if a cell satisfied this condition.')
        end
    end
end
metric = cell2mat(Vals);
[conditions, conditions_id] = findgroups(Condition_id); % remove empty conditions
days_elapsed = cell2mat(Days_elapsed);
n_electrodes = cell2mat(N_electrodes);
%% normalize by electrodes
if P_in.normalize_by_electrodes
    metric=metric./n_electrodes;
    y_label = [y_label '/electrode'];
else
    if ~strcmp(P_in.metric, 'Vpp')
        warning('Because of conditioning, it is recommended to normalize the number of electrodes')
    end
end
%% Plot
figure('Position', P.figure_position_longevity);
set(gca, P.axes_properties_behavior{:})
clear boots
for i = 1:numel(unique(conditions))
    boots = nan(n_boots,numel(t_bin_centers));
    for j=1:numel(t_bin_centers)
        idx = days_elapsed>t_bin_edges(j) & days_elapsed<t_bin_edges(j+1) & ...
              conditions==i;
        if sum(idx) > 1
            boots(:,j) = bootstrp(n_boots,@mean,metric(idx(:)));
        elseif sum(idx)==1
            boots(1:n_boots,j) = metric(idx(:));
        else
            continue
        end
    end
    hdl(i)=shadedErrorBar(t_bin_centers,boots,{@mean,@std});hold on;
    hdl(i).mainLine.LineWidth=1;
    hdl(i).mainLine.Marker = 'o';
    hdl(i).patch.FaceAlpha=P_in.FaceAlpha;
    if numel(unique(conditions)) == 1
        the_color = [0,0,0];
    else
        the_color = P.color_order(i,:);
    end
    hdl(i).mainLine.Color = the_color;
    hdl(i).patch.FaceColor = the_color;
end
if numel(unique(conditions))>1
    legend_labels = arrayfun(@(x) [char(string(x))], conditions_id, 'uni', 0);
    legend([hdl.mainLine],legend_labels);
end
set(gca,'xscale','log',...
        'xlim',[min(t_bin_edges), max(t_bin_edges)],...
        'xtick',t_bin_centers);
ylim([0,1].*ylim)
xlabel('Days since implant')
ylabel(y_label); 
%% Save
for i = 1:numel(P.figure_image_format)
    file_name = ['average_' P_in.metric];
    if ~isempty(P_in.condition_on)
        file_name = [file_name '_conditioning_on_' P_in.condition_on];
    end
    saveas(gcf, [P.plots_folder filesep file_name], P.figure_image_format{i})
end