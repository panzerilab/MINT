%%% *function R_1d = map_Nd_array_to_1d(R_Nd)*
%%%
%%% ### Description
%%% This function maps a N-dimensional array defined by `R_Nd` to 1D vector.
%%% The function is useful when analysing multi-dimensional activity (either
%%% in time on a single unit or across units) when interested in understanding
%%% if specific activation patterns result in different information levels.
%%%
%%% The multi-dimensional input is mapped on a single dimensional vector of
%%% the same size as the possible combinations of all responses in the
%%% input multi-dimensional activity.
%%%
%%% As an example, the 2-dimensional input `R_nd` below:
%%%
%%%     R_nd = [ 1 2 1; 2 3 1]
%%%
%%% Is mapped on the rows of the matrix of all possible 2-dimensional responses:
%%%
%%%     R_nd = [ 1 1; 2 1; 1 2; 2 2; 1 3; 2 3]
%%%
%%% The output `R_1d` is then:
%%%
%%% 	R_nd = [3 6 1]
%%%
%%% ### Inputs:
%%% - *R_Nd*: *nDims X nTrials* input response matrix
%%%
%%% ### Outputs:
%%% - *R_1D*: output 1D response based on the N-dimensional space defined by `R_Nd`
%%%
