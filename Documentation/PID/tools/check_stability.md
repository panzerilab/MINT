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
