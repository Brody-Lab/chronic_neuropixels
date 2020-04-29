function [theSEM] = sem(x, varargin)


theSEM = nanstd(x,varargin{:}) ./ sqrt(sum(~isnan(x)));