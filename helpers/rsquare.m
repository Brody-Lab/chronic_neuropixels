function R2 = rsquare(y,yhat)
    % PURPOSE:  calculate r square using data y and estimates yhat
    % -------------------------------------------------------------------
    % USAGE: R2 = rsquare(y,yhat)
    % where: 
    %        y are the original values as vector or ND arrays and
    %        yhat are the estimates calculated from y using a regression, given in
    %        the same form (vector or raster) as y
    % -------------------------------------------------------------------------
    % OUTPUTS:
    %        R2 is the r square value calculated using 1-SS_E/SS_T (by
    %        default, with other options)
    % -------------------------------------------------------------------
    % Note: NaNs in either y or yhat are deleted from both sets.
    %
    % Felix Hebeler, Geography Dept., University Zurich, Feb 2007
    %
    % can calculate R2 using a variety of methods
    %
    % edited by AGB, 2015
    
    verbose=false;    
    method=1;
    %% check args
    if nargin ~= 2
        error('This function needs exactly 2 input arguments!');
    end
    sz1=size(y);
    sz2=size(yhat);
    if any(sz1~=sz2)
        error('Inputs have different sizes.');
    end
    if verbose && numel(y)==1
        mssg(0,'Inputs are scalars!','color',[1 0 0]);
    end
    %% reshape
    yhat=yhat(:);
    y=y(:);
    %% delete NaNs
    while any(isnan(y)) || any(isnan(yhat))
        if sum(isnan(y)) >= sum(isnan(yhat)) 
            yhat(isnan(y))=[];
            y(isnan(y))=[];
        else
            y(isnan(yhat))=[]; 
            yhat(isnan(yhat))=[];
        end
    end
    switch method
        case 1 % 1 - SSe/SSt (DEFAULT)
            R2 = 1 - ( sum( (y-yhat).^2 ) / sum( (y-mean(y)).^2 ) );
        case 2 % SSr/SSt
            R2 = sum((yhat-mean(y)).^2) / sum( (y-mean(y)).^2 ) ;
        case 3 % squared Pearson correlation
            R2 = corr(yhat,y).^2;
        otherwise
            error('Unknown method');
    end
    %%
    if verbose && ( R2<0 || R2>1 )
        mssg(0,['R^2 of ',num2str(R2),' : yhat does not appear to be the estimate of y from a within-fold regression.'])
    end
end