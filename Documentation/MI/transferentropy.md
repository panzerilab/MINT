%%% *function transferentropies = transferentropy(X, Y, opts, outputsList)*
%%%
%%% ### Description
%%% Computes transfer entropy from signal $X$ to signal $Y$.
%%%
%%% ### Inputs:
%%% - *X*: *nTimePts X nTrials* matrix containing the "causing" signal data matrix.
%%% - *Y*: *nTimePts X nTrials* matrix containing the "caused" signal data matrix.
%%% - *opts*: options structure.
%%% - *outputsList*: cell array of char arrays of strings specifying the quantities that the function is computing. These are going to be returned, in the same order as specified in this list in *entropies*. See further notes section for more details on the available options.
%%%
%%% ### Outputs:
%%% - *transferentropies*: cell array of same length as *outputsList* returning the specified outputs.
%%%
%%% ### Further notes
%%% ### The data matrices
%%%
%%% The function requires that the user specifies two types of data: the
%%% "causing" signal X and the "caused" signal Y. This toolbox allows to
%%% estimate transfer entropy, i.e. how much the knowledge of the past of X
%%% allows to increase the predictability over the present of Y.
%%% 
%%% Both X and Y must be matrices of size *nTimePts X nTrials* where **nTimePts* is the
%%% number of points recorded in each trial while *nTrials* is the number of
%%% trials. Trials can be different recordings of the two signals or they
%%% can just reflect the way the user wishes to break the data into blocks
%%% for statistical testing (see [Bootstrapping](#bootstrapping) section.
%%% 
%%% It is important to remind that, in order to compute transfer entropy
%%% trials of X and corresponding trials of Y must be simultaneously
%%% recorded.
%%%
%%% #### The options structure
%%% The function `transferentropy` is essentially a wrapper around function
%%% [`entropy`](MI/entropy).
%%% Consequently, `transferentropy` options include many of the options to
%%% `entropy` (exceptions are represented by the number-of-trials per
%%% stimulus opts.nt array and by the usage of the bootstrap option opts.btsp
%%% option, see [Bootstrapping](#bootstrapping) section.
%%%
%%% Below is a list of the required and optional option parameters: for a 
%%% description of the inherited input options the user can refer to the
%%% [documentation](MI/BinnedMethods/entropy) of function `entropy`.
%%%
%%% ##### `transferentropy`-specific options
%%% - opts.taux (*required*): lags considered for the causing signal. Lags are specified as arrays of **strictly negative** integers.
%%% - opts.tuay (*required*): lags considered for the caused signal. Lags are specified as arrays of **strictly negative** integers.
%%% - opts.trperm (*optional*): number of random trial permutations. **Cannot be specified** in conjunction with *opts.btsp* option.
%%%
%%% ##### Options inherited from function `entropy`
%%% - opts.method (*required*): probability estimation method.
%%% - opts.bias (*required*): bias correction method.
%%% - opts.btsp (*optional*): number of bootstrap iterations. **Cannot be specified** in conjunction with *opts.trperm* option.
%%% - opts.verbose (*optional*): perform extensive checks on inputs, default true.
%%%
%%% #### The output list
%%% To specify which IT quantities need to compute, one or more of the following strings has to be specified:
%%%
%%% | Option  | Description                                                   |
%%% |---------|---------------------------------------------------------------|
%%% | 'TE'    | $X \to Y$ transfer entropy                                    |
%%% | 'TEsh'  | $X \to Y$ transfer entropy with shuffle correction            |
%%% | 'NTE'   | normalized $X \to Y$ transfer entropy                         |
%%% | 'NTEsh' | normalized $X \to Y$ transfer entropy with shuffle correction |
%%%
%%%  Outputs are returned IN THE SAME ORDER as that specified in the output list.
%%%
%%% #### Trial-permutation and bootstrapped estimates
%%% For `opts.trperm > 0`, *trial-permutation* estimates are computed after 
%%% trials in X and Y have been randomly permuted. This can be useful for 
%%% bias performing statistics on the transfer-entropy values or to estimate
%%% the effect of a causing signal which is common to both X and Y. During
%%% trial-permutation the time-consistency of *X* and *Y* is kept, while
%%% trials are shuffled.
%%%
%%% *Bootstrapped estimates*, instead, destroy the time-consistency within
%%% trials by randomly shuffling time points within the same time lag
%%% considered. Bootstrap estimates are returned if `opts.btsp > 0`.
%%% 
%%% Trial-permutation and bootstrapping are *mutually exclusive*.
%%% Trial-permuted and bootstrap estimates can be returned for all
%%% [outputs](the-output-list) of the function.
%%%
%%% Bootstrap (or trial-permuted) estimates are returned to the user in
%%% the form of an array of length `opts.trperm` (or `opts.btsp`) concatenated
%%% to the actual transfer-entropy estimate. For example, suppose that
%%% TE is computed with Opt.btsp = 20. In this case the output
%%% corresponding to TE will be an array of length 21: the first
%%% element is the actual entropy estimate while the remaining 20
%%% elements correspond to 20 distinct bootstrap estimates.
%%%
%%%
%%%                Actual | bootstrap (trial-permuted)
%%%              Estimate | estimates
%%%                       |
%%%           Index:    1 | 2   3   4      Opt.btsp+1 = 21
%%%                   ----|------------- ... -----
%%%                   | x | x | x | x |      | x |
%%%                   ----|------------- ... -----
%%%                       |
%%%
