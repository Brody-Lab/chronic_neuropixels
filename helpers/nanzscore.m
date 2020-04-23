function [z,mu,sigma] = nanzscore(x,flag,dim,varargin)
%ZSCORE Standardized z score.
%   Z = NANZSCORE(X) returns a centered, scaled version of X, the same size as X.
%   For vector input X, Z is the vector of z-scores (X-NANMEAN(X)) ./ NANSTD(X). For
%   matrix X, z-scores are computed using the mean and standard deviation
%   along each column of X.  For higher-dimensional arrays, z-scores are
%   computed using the mean and standard deviation along the first
%   non-singleton dimension. NaNs are treated as missing values.
%
%   The columns of Z have sample mean zero and sample standard deviation one
%   (unless a column of X is constant, in which case that column of Z is
%   constant at 0).
%
%   [Z,MU,SIGMA] = NANZSCORE(X) also returns NANMEAN(X) in MU and NANSTD(X) in SIGMA.
%
%   [...] = NANZSCORE(X,1) normalizes X using NANSTD(X,1), i.e., by computing the
%   standard deviation(s) using N rather than N-1, where N is the length of
%   the dimension along which ZSCORE works.  NANZSCORE(X,0) is the same as
%   NANZSCORE(X).
%
%   [...] = NANZSCORE(X,FLAG,DIM) standardizes X by working along the dimension
%   DIM of X. Pass in FLAG==0 to use the default normalization by N-1, or 1
%   to use N.
%
%   [...] = NANZSCORE(X,...,'balance',groups) applies "balanced" z-scoring
%   following the method proposed in Kang & Maunsell (2012). J.
%   Neurophysiol. for grouped data, where groups is a logical array specifying group IDs.
%   (This is useful for calculating grand CP across stimulus conditions with
%   non-uniform choice distributions.)

%   Adrian Bondy, 2015
%   based on Matlab z-score, copyright 1993-2006 The MathWorks, Inc. 

% [] is a special case for std and mean, just handle it out here.
if isempty(x), z = []; return; end

if nargin < 2
    flag = 1;
end
if nargin < 3
    % Figure out which dimension to work along.
    dim = find(size(x) ~= 1, 1);
    if isempty(dim), dim = 1; end
end

% Compute X's mean and sd, and standardize it
groups=[];
balanceInds = strncmp(varargin,'balance',4);
if ~isempty(balanceInds)
    groups=varargin{find(balanceInds)+1};
    if ~islogical(groups) || any(size(groups)~=size(x))
        error('Groups must be a logical array of the same size as x.');
    end
    counts = Counts(groups);
    if numel(counts)==1 || min(counts)<2
        warning('nanzscore:GroupsTooUnevenForBalanced','Will use unbalanced z-scoring because there are only enough samples (at least 2) from a single group.');
        groups=[];
    end    
end
stdFun = @(x)nanstd(x,flag,dim);
meanFun = @(x)nanmean(x,dim);
if isempty(groups)
    mu = meanFun(x);
    sigma = stdFun(x);    
else
    group1vals = x(groups);
    group2vals = x(~groups);
    m1 = meanFun(group1vals);
    m2 = meanFun(group2vals);
    mu = ( m1 + m2 ) / 2;    
    sigma = sqrt ( ( stdFun(group1vals).^2 + stdFun(group2vals).^2  + (m1-m2).^2/2 ) / 2 ) ;
end
sigma0 = sigma;
sigma0(sigma0==0) = 1;
z = bsxfun(@minus,x, mu);
z = bsxfun(@rdivide, z, sigma0);

