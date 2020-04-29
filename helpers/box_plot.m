function box_plot(x,y,varargin)
    p=inputParser;
    p.addRequired('x',@(x)validateattributes(x,{'numeric'},{'vector'}));
    p.addRequired('y',@(d)validateattributes(d,{'numeric'},{'2d','ncols',numel(x),}));
    p.addParameter('labels',@(x)validateattributes(x,{'char'},{'size',size(x)}));
    p.addParameter('ylim',[],@(x)validateattributes(x,{'numeric'},{'numel',2,'increasing'}));
    p.addParameter('center',@median,@(x)validateattributes(x,{'function_handle'},{'scalar'}));
    p.addParameter('whisker_prctile_range',[2.5 97.5],...
        @(x)validateattributes(x,{'numeric'},{'numel',2,'increasing','nonnegative','<=',100}));
    p.addParameter('box_prctile_range',[25 75],...
        @(x)validateattributes(x,{'numeric'},{'numel',2,'increasing','nonnegative','<=',100}));
    p.addParameter('axes',gca,@(x)validateattributes(x,{'matlab.graphics.axis.Axes'},{'nonempty'}));
    p.addParameter('boxWidth',0.15,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('boxEdgeWidth',1,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('whiskerWidth',2,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('centerWidth',2,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('centerLength',0.2,@(x)validateattributes(x,{'numeric'},{'positive','scalar'}));
    p.addParameter('boxEdgeColor',[1 1 1]/2,@(x)validateattributes(x,{'numeric'},{'vector','numel',3}));
    p.addParameter('boxFaceColor',[1 1 1]/2,@(x)validateattributes(x,{'numeric'},{'vector','numel',3}));
    p.addParameter('centerColor',[0 0 0],@(x)validateattributes(x,{'numeric'},{'vector','numel',3}));
    p.addParameter('whiskerColor',[0 0 0],@(x)validateattributes(x,{'numeric'},{'vector','numel',3}));
    p.parse(x,y,varargin{:});
    params=p.Results;
    axes(params.axes);
    n_x = numel(x);
    for i=1:n_x
       patch([-1 1 1 -1 -1]*params.boxWidth + x(i),...
           prctile(y(:,i),params.box_prctile_range([1 1 2 2 1])),...
           params.boxFaceColor,'linewidth',params.boxEdgeWidth,'FaceColor',params.boxFaceColor,...
           'EdgeColor',params.boxEdgeColor);
       line([1 1]*x(i),prctile(y(:,i),params.whisker_prctile_range),'linewidth',params.whiskerWidth,...
           'color',params.whiskerColor);
       line([-1 1]*params.centerLength+x(i),[1 1]*params.center(y(:,i)),'linewidth',params.centerWidth,...
           'color',params.centerColor);
    end
    set(gca,'xtick',x,'xticklabel',params.labels);
end