%%% *function [X, nt] = buildx(Y, X_in)*
%%%
%%% ### Description
%%% Given input 2D matrix *X_in* (typically neural response) and 1D input vector *Y* (typically stimuli array), it reshapes the X matrix to a 3D shape *nDimsX X max(nTrials) X nY*.
%%%
%%% ### Inputs:
%%% - *Y*: *1 X nTrials* 1D array.
%%% - *X_in*: *nDimsX X nTrials* 2D matrix describing the response of each of the *nDimsX* dimensions for each trial.
%%%
%%% ### Outputs:
%%% - *X*: *nDims X max(nTrials) X nY* reshaped input *X_in* matrix, where *max(nTrials)* is the maximum across all the number of trials for each discrete value of Y.
%%% - *nt*: *nY X 1* array reporting the number of trials for each of the discrete values of *Y*.
%%%
