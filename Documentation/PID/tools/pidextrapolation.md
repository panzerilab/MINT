%%% function [PID_corrected, PID_naive] = pidextrapolation(Y_b, X1_b, X2_b, opts, btspflag)
%%%
%%% ### Description
%%% pidextrapolation computes the partial information decomposition (PID) for data variables, corrected for estimation biases that might occur in smaller sample sizes. 
%%% The function employs bootstrap and extrapolation methods to estimate the bias and correct the initial PID values.
%%%
%%% ### Inputs:
%%% - Y_b, X1_b, X2_b: Data vectors representing the target variable (Y_b) and two source variables (X1_b, X2_b), respectively.
%%% - opts: A structure containing options for the PID computation, including:
%%%     - xtrp: The number of extrapolation points.
%%%     - n_trials: Number of trials or samples in the data.
%%%     - function: Function handle for computing PID.
%%%     - parallel: Flag to indicate if parallel processing should be used.
%%%     - bias: Type of bias correction to apply ('qe' for quadratic extrapolation, 'le' for linear extrapolation).
%%%     - btspflag: A flag indicating whether bootstrap resampling should be performed.
%%%
%%% ### Outputs:
%%% - PID_corrected: A vector containing the bias-corrected PID values for the data.
%%% - PID_naive: A vector containing the initially computed PID values without correction for bias.
