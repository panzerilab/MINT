function [PID_v1, PID_v21, PID_v22] = check_stability(Y, X1, X2, opts)
%%% *function [PID_v1, PID_v21, PID_v22] = check_stability(Y, X1, X2, opts)*
%%%
%%% ### Description
%%% This function assesses the stability of Partial Information Decomposition (PID)
%%% results by calculating the PID on the full dataset and two equal partitions of it. 
%%% This helps to understand the stability of PID measurements across different dataset splits.
%%%
%%% ### Inputs:
%%% - *Y*: must be an array of *nDimensions X nTrials* elements representing the first signal
%%% - *X1*: must be an array of *nDimensions X nTrials* elements representing the second signal
%%% - *X2*: must be an array of *nDimensions X nTrials* elements representing the third signal
%%% - *opts*: options used to calculate PID (see further notes).
%%%
%%% ### Outputs:
%%% - *PID_v1*: PID results for the entire dataset.
%%% - *PID_v21*: PID results for the first partition of the dataset.
%%% - *PID_v22*: PID results for the second partition of the dataset.
%%%   Each output contains fields:
%%%   - `qe`: bias corrected estimate of the PID.
%%%   - `naive`: Naive estimate of the PID.
%%%
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | field and default value             | description                                                                                           | Possible values
%%% |-------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%% | opts.bin_methodX1 = `'eqpop'`       | binning method for the X1 signal                                                                      | 'none', 'eqpop', 'eqspace', 'ceqspace', 'gseqspace'                                                                                                          |
%%% | opts.bin_methodX2 = `'eqpop'`       | binning method for the X2 signal                                                                      | 'none', 'eqpop', 'eqspace', 'ceqspace', 'gseqspace'                                                                                                          |
%%% | opts.bin_methodY = `'none'`         | binning method for the Y signal                                                                       | 'none', 'eqpop', 'eqspace', 'ceqspace', 'gseqspace'                                                                                                          |
%%% | opts.n_binsX1 = 3                   | number of bins to reduce X1 dimensionality                                                            | int > 1                                                                                                                                                      |
%%% | opts.n_binsX2 = 3                   | number of bins to reduce X2 dimensionality                                                            | int > 1                                                                                                                                                      |
%%% | opts.n_binsY = N/A                  | number of bins to reduce Y  dimensionality                                                            | int > 1                                                                                                                                                      |
%%% | opts.btsp = 0                       | number of bootstrap operations to perform for significance testing                                    | int >= 0                                                                                                                                                     |
%%% | opts.btsp_variables = N/A           | list of variables to be bootstrapped, specified as a cell array of strings.                           | cell array containing one (or more) strings (`"X1"`, `"X2"` or `"Y"`)                                                                                        |
%%% | opts.btsp_type = N/A                | type of bootstrapping to be applied for each variable                                                 | cell array of same size of opts.btsp_variables, each element contains the shuffling type (either 'all', 'X1conditioned', 'X2conditioned', or 'Yconditioned') |
%%% | opts.old_output = 0                 | flag to indicate if the output will follow the old output from the NIT toolbox (1 if yes)             | int 0 or 1                                                                                                                                                   |
%%% | opts.xtrp = 50                      | how many iterations repetitions of the extrapolation procedure should be performed.                   | int > 0                                                                                                                                                      |
%%% | opts.redundancy_measure = 'I_BROJA' | name of the measure of the redundancy between sources                                                 | 'I_BROJA', 'I_min' or 'I_MMI'                                                                                                                                |
%%% | opts.function = @pidBROJA           | function handle of the redundancy measure between sources, not necessary if redundancy_measure is set | any function handle that accepts a probability distribution as input                                                                                         |
%%% | opts.parallel = 0                   | indicates if the extrapolation and bootstrap procedures run in parallel (1 if yes)                    | int 0 or 1                                                                                                                                                   |
%%%
%%% Binning methods
%%% - 'none' (no binning)
%%% - 'eqpop' (evenly populated binning)
%%% - 'eqspace' (evenly spaced binning)
%%% - 'ceqspace' (centered evenly spaced binning)
%%% - 'gseqspace' (gaussian evenly spaced binning)




% X1X2cS = buildx(Y, X1, X2);
n_trials = consistency_check({Y X1 X2});
allowedMethods = ["none", "eqpop", "eqspace", "ceqspace", "gseqspace"];
binMethods = {opts.bin_methodX1, opts.bin_methodX2, opts.bin_methodY};
methodNames = {'bin_methodX1', 'bin_methodX2', 'bin_methodY'};
for i = 1:length(binMethods)
    currentMethod = binMethods{i};
    currentName = methodNames{i};
    if ~any(arrayfun(@(x) strcmp(x, currentMethod), allowedMethods))
        if ~exist(currentMethod, 'file')
            msg = "User defined binning method: `" + currentMethod + "` for " + currentName + " cannot be found in the current MATLAB search path.";
            error('PID:InvalidBinningMethod', msg);
        elseif exist(currentMethod, 'file') ~= 2
            msg = "User defined binning method: `" + currentMethod + "` for " + currentName + " could be found in the current MATLAB search path but does not appear to be a MATLAB function.";
            error('PID:InvalidBinningMethod', msg);
        end
    end
end
%opts = check_opts(opts);
opts.n_trials =n_trials;
opts.bias = 'qe';

if ~strcmp(opts.bin_methodX1, 'none')
    X1_b = binr(X1, opts.n_binsX1, opts.bin_methodX1)+1;
else
    X1_b = X1+1;
end
if ~strcmp(opts.bin_methodX2, 'none')
    X2_b = binr(X2, opts.n_binsX2, opts.bin_methodX2)+1;
else
    X2_b = X2+1;
end
if ~strcmp(opts.bin_methodY, 'none')
    Y_b  = binr(Y, opts.n_binsY, opts.bin_methodY)+1;
else
    Y_b  = Y+1;
end
if length(X1_b(:,1)) > 1
    X1_b = reduce_dim(X1_b, 1);
end
if length(X2_b(:,1)) > 1
    X2_b = reduce_dim(X2_b, 1);
end
if length(Y_b(:,1)) > 1
    Y_b = reduce_dim(Y_b, 1);
end

ri = randperm(n_trials, n_trials);
Y_bs = Y_b(ri);
X1_bs = X1_b(ri);
X2_bs = X2_b(ri);

[uniqStim,ia,ic] = unique(Y_bs);

partition1 = [];
partition2 = [];
for stim=1:length(uniqStim)
    mask= (Y_bs==uniqStim(stim));
    Ytotstim = Y_bs(mask);
    X1totstim = X1_bs(mask);
    X2totstim = X2_bs(mask);

    bin_edge = round(length(Ytotstim)/2);
    partition1 = [partition1 [Ytotstim(1,1:bin_edge); X1totstim(1,1:bin_edge); X2totstim(1,1:bin_edge)]];
    partition2 = [partition2 [Ytotstim(1,bin_edge+1:end); X1totstim(1,bin_edge+1:end); X2totstim(1,bin_edge+1:end)]];
end

opts.bin_methodY  = 'none';
opts.bin_methodX1 = 'none';
opts.bin_methodX2 = 'none';

[PID_v1qe, PID_v1naive, ~, ~, ~, ~]= PID(Y_b, X1_b, X2_b, opts);
[PID_v21qe, PID_v21naive, ~, ~, ~, ~]= PID(partition1(1,:), partition1(2,:), partition1(3,:), opts);
[PID_v22qe, PID_v22naive, ~, ~, ~, ~]= PID(partition2(1,:), partition2(2,:), partition2(3,:), opts);


PID_v1.qe = PID_v1qe;
PID_v1.naive = PID_v1naive;

PID_v21.qe = PID_v21qe;
PID_v21.naive = PID_v21naive;

PID_v22.qe = PID_v22qe;
PID_v22.naive = PID_v22naive;
end
