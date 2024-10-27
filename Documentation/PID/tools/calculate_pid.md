%%% Given three variables Y, X1, and X2, calculate_pid calculates the Partial Information Decomposition (PID)
%%% of X1 and X2 about Y.
%%%
%%% Inputs:
%%% - Y: Array containing values for variable Y
%%% - X1: Array containing values for variable X1
%%% - X2: Array containing values for variable X2
%%% - n_trials: Total number of trials in the dataset
%%% - n_split_trials: Number of trials to be randomly sampled for each calculation
%%%
%%% Output
%%% - PID value for X1 and X2 about Y
%%%
%%%
%%% random sample n_split_trials from data in each of the
%%% n_draws and calculate II on the reduced dataset
%%% random sample n_split_trials from data in each of the
%%% n_draws and calculate II on the reduced dataset
