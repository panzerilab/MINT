%%% *function binnedX = eqpop(X, nb, varargin)*
%%%
%%% ### Description
%%% Builds edges for equipopulated binning. Unlike the other binning functions provided, that bin data according to the values of X, equally-populated binning bins according to the count. For this reason, this function returns already the binned version of X and **not** the bin edges.
%%%
%%% ### Inputs:
%%% - *X*: array of values to be discretized.
%%% - *nb*: number of bins.
%%%
%%% ### Outputs:
%%% - *binnedX*: unlike the other binning functions, that bin.
