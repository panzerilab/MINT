function PID = calculate_pid(Y,X1,X2,n_trials,n_split_trials)
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

if nargin < 5
    msg = 'not enough input arguments.';
    error('calculatepid:notEnoughInput', msg);
end

variables = {{Y, 'Y'}, {X1, 'X1'}, {X2, 'X2'}, {n_trials, 'n_trials'}, {n_split_trials, 'n_split_trials'}};
nanVariableName = '';
for i = 1:length(variables)
    var = variables{i}{1};
    name = variables{i}{2};
    if any(isnan(var), 'all')
        nanVariableName = name;
        break;
    end
end

if ~isempty(nanVariableName)
    msg = sprintf('%s contains NaNs. Aborting.', nanVariableName);
    error('calculatepid:NaNInput', msg);
end


ri = randperm(n_trials, n_split_trials);
Yr = Y(ri);
X1r = X1(ri);
X2r = X2(ri);
[prob_matrix, ~, nY, nX1, nX2, n_s_d] = build_p(Yr', X1r', X2r');
PID = pidBROJA(prob_matrix);
end