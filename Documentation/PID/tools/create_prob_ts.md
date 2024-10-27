%%% function p_ts = create_prob_ts(p_distr, sources)
%%%
%%% ### Description
%%% Function to create a bivariate prob distribution between target and sources. If sources is just one variable, it's the same as marginalizing
%%% the other variables. If, for example, the target is X with size 2 and the sources are two variables Y and Z with number of values 3 and 4, the
%%% resulting matrix is (|X|, |Y|*|Z|)
%%%
%%% ### Inputs:
%%% - p_distr: A multidimensional array representing the original probability distribution.
%%% - sources: A vector indicating the indices of the dimensions in p_distr that should be treated as sources.
%%%
%%% ### Outputs:
%%% - p_ts: A two-dimensional matrix representing the joint probability distribution between
%%%         the target and source variables.
