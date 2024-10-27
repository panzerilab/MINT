%%% *function entropies = entropy(X, Y, opts, outputsList)*
%%%
%%% ### Description
%%% Compute entropy and entropy-like quantities.
%%%
%%% ### Inputs:
%%% - *X*: *nDimensionsX X nTrials* array including input data (typically the neural response).
%%% - *Y*: *nDimensionsY X nTrials* array including input data (typically the stimulus). In principle the values of Y can also represent a multi-dimensional neural response (if the MI between two neuronal populations is of interest). While the MI calculated is symmetrical to swapping X and Y, the user should be aware that the information breakdown quantities are calculated only on X. 
%%% - *opts*: options structure (see further notes section for more details).
%%% - *outputsList*: cell array of char arrays of strings specifying the quantities that the function is computing. These are going to be returned, in the same order as specified in this list in *entropies*. See further notes section for more details on the available options.
%%%
%%% ### Outputs:
%%% - *entropies*: cell array of same length as *outputsList* returning the specified outputs.
%%%
%%% ### Further notes
%%% #### The options structure
%%% The options structure can include any the following fields:
%%%
%%% - opts.n_binsX: number of bins array to be used on X (can be a scalar or a vector of `size(nb) = size(X,1)`). If a scalar is specified then all dimensions are binned using the same number of bins. Alternatively, each dimension is binned using its specific number of bins. 
%%% - opts.n_binsY: number of bins array to be used on Y (can be a scalar or a vector of `size(nb) = size(X,1)`). If a scalar is specified then all dimensions are binned using the same number of bins. Alternatively, each dimension is binned using its specific number of bins. 
%%% - opts.bin_methodX: binning method for X (can be a scalar or a cell array of `size(method) = size(X,1)`). If a scalar is specified then all dimensions are binned using the same method. Alternatively, each dimension is binned using its specific method. If not specified it is assumed that no binning is necessary and X is already discrete.
%%% - opts.bin_methodY: binning method for Y (can be a scalar or a cell array of `size(method) = size(Y,1)`). If a scalar is specified then all dimensions are binned using the same method. Alternatively, each dimension is binned using its specific method. If not specified it is assumed that no binning is necessary and Y is already discrete.
%%% 
%%% Possible options for binning methods are:
%%%
%%% | Option       | Description                        |
%%% |--------------|------------------------------------|
%%% | `'none'`     | No binning                         |
%%% | `'eqpop'`    | Equally populated binning          |
%%% | `'eqspace'`  | Equispaced binning                 |
%%% | `'ceqspace'` | Centered equispaced spaced binning |
%%% | `'geqspace'` | Gaussian equispaced spaced binning |
%%%
%%% see the documentation of [`binr`](tools/Binning/binr) function for more details on each of these binning strategies.
%%%
%%% - opts.method: this field specifies which estimation method to use and can be one of the following strings:
%%%
%%% | Option | Description     |
%%% |--------|-----------------|
%%% | `'dr'` | Direct method   |
%%% | `'gs'` | Gaussian method |
%%%
%%% The direct method requires X values to be discretized into non-negative integer values. This is handled internally to the code through the *opts.bin_methodX* option. It can also be done by the user prior to the call to this function but the user should be aware that the binned values should be strincly non-negative. **Failing to properly discretizing the X array will result in Matlab crashing.**
%%%
%%% - opts.bias: this field specifies the bias correction procedure. It can be one of the following strings:
%%%
%%% | Option    | Description             |
%%% |-----------|-------------------------|
%%% | `'qe'`    | Quadratic extrapolation |
%%% | `'pt'`    | Panzeri & Treves 1996   |
%%% | `'bub'`   | Paninski 2003           |
%%% | `'gsb'`   | Gaussian bias           |
%%% | `'naive'` | Biased naive estimates  |
%%%
%%% - opts.btsp (optional, default: *opt.btsp = 0*): this field must be a (non-negative) scalar specifying how many bootstrap estimates to compute.These estimates are useful as an additional bias correction or for performing statistics on the entropy values.
%%%       Bootstrap estimates are performed by means of pairing X and Y at random and computing the entropy quantities for these random pairings; each estimate corresponds to a different random pairing configuration.
%%%       See the examples below for additional information on how to use this option.
%%% - opt.xtrp (optional, default: *opt.xtrp = 0*): this field must be a (non-negative) scalar specifying how many  iterations repetitions of the extrapolation procedure should be performed. Extrapolations procedure (such as the quadratic extrapolation) perform bias estimation by computing entropy values on sub-groups of the available trials. These subgroups are created randomly. The xtrp option allows to average the extrapolation values over as many different random partitions as specified by the parameter.
%%% - opts.verbose (optional, default *opt.verbose = true*): if this field exists and is set to true a summary of the selected options is displayed and additional checks are performed on the input variables. No warnings are displayed unless this options is enabled. This feature is useful to check whether INFORMATION is being called correctly. It is therefore highly reccomended for new users or when first running of the program with new input options. However, keep in mind that these checks drammatically increases computation time and are thus not reccommended for computationally intensive session. If a custom bias correction function is called (see "BUILDING AND CALLING CUSTOM BIAS CORRECTION FUNCTIONS" below) tests are performed on the invoked function to check whether it satisfied the requirements of the toolbox.
%%%   
%%% #### The output list
%%% To specify which IT quantities need to compute, one or more of the following strings has to be specified:
%%%
%%% | Option  | Description       |
%%% |---------|-------------------|
%%% | 'HX'    | $H(X)$            |
%%% | 'HXY'   | $H(X|Y)$          |
%%% | 'HlX'   | $H_{lin}(X)$      |
%%% | 'HiX'   | $H_{ind}(X)$      |
%%% | 'HiXY'  | $H_{ind}(X|Y)$    |
%%% | 'ChiX'  | $Chi(X)$          |
%%% | 'HshX'  | $H_{sh}(X)$       |
%%% | 'HshXY' | $H_{sh}(X|Y)$     |
%%%
%%%  Outputs are returned IN THE SAME ORDER as that specified in the output list.
%%% IMPORTANT: Not all combinations of method, bias and output options are possible. For example, bias correction 'pt' can only be used together  with method 'dr'. The allowed combinations of method, bias and output options are summarized in the following tables: 
%%%   
%%% #### Allowed outputs/bias combinations for direct method estimation
%%% Legend:
%%% - **X**: combination available
%%% - **n**: naive estimate returned
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'bub' | 'gsb' | 'user defined' |
%%% |---------|---------|-------|-------|-------|-------|-----------------
%%% | 'HX'    |    X    |   X   |   X   |   X   |   n   |       X        |
%%% | 'HXY'   |    X    |   X   |   X   |   X   |   n   |       X        |
%%% | 'HlX'   |    X    |   X   |   X   |   X   |   n   |       X        |
%%% | 'HiX'   |    X    |   X   |   n   |   n   |   n   |       n        |
%%% | 'HiXY'  |    X    |   X   |   X   |   n   |   n   |       n        |
%%% | 'ChiX'  |    X    |   X   |   n   |   n   |   n   |       n        |
%%% | 'HshX'  |    X    |   X   |   X   |   X   |   n   |       X        |
%%% | 'HshXY' |    X    |   X   |   X   |   X   |   n   |       X        |
%%% 
%%% #### Allowed outputs/bias combinations for Gaussian method estimation
%%% Legend:  
%%% - **X**: combination available
%%% - **n**: naive estimate returned
%%% - **0** : zero returned
%%%
%%% |         | 'naive' | 'qe'  | 'pt'  | 'bub' | 'gsb' | 'user defined' |
%%% |---------|---------|-------|-------|-------|-------|-----------------
%%% | 'HX'    |    X    |   X   |   n   |   n   |   X   |       n        |
%%% | 'HXY'   |    X    |   X   |   n   |   n   |   X   |       n        |
%%% | 'HlX'   |    X    |   X   |   n   |   n   |   X   |       n        |
%%% | 'HiX'   |    0    |   0   |   0   |   n   |   0   |       n        |
%%% | 'HiXY'  |    X    |   X   |   n   |   n   |   X   |       n        |
%%% | 'ChiX'  |    0    |   0   |   0   |   n   |   0   |       n        |
%%% | 'HshX'  |    X    |   X   |   n   |   n   |   X   |       n        |
%%% | 'HshXY' |    X    |   X   |   n   |   n   |   X   |       n        |
%%%
%%% #### Output options with bootstrap
%%% Bootstrap estimates make sense (and can thus be computed) only for
%%% the following output quantities: 'HXY', 'HiX', 'HiXY', 'ChiX' and
%%% HshXs'.
%%%
%%% If the user only specifies the number of bootstrap repetitions
%%% (through the BTSP parameter) bootstrap estimates are computed for
%%% any of the above quantities which appears in the output list.
%%%
%%% However users are given the opportunit to precisely select which
%%% quantities to compute bootstrap for by means of appending 'bs' to
%%% the output string name of an output quantitie as follows:
%%%
%%% | Quantity | Bootstrapped quantity  |
%%% |----------|------------------------|
%%% | 'HiX'    | 'HiXbs'                |
%%% | 'HiXY'   | 'HiXYbs'               |
%%% | 'ChiX'   | 'ChiXbs'               |
%%% | 'HshXY'  | 'HshXYbs'              |
%%%
%%%  Carefully selecting for which quantities to compute bootstrap can
%%%  greately decrease computation times. For example, bootstrap
%%%  estimates of HiX and ChiX are often useless although their
%%%  computation can be very time consuming.
%%%
%%%  Bootstrap estimates are returned to the user in the form of an
%%%  array of length Opt.btsp concatenated to the actual entropy
%%%  estimate. For example, suppose that HXY is computed with Opt.btsp =
%%%  20. In this case the output corresponding to HXY will be an array
%%%  of length 21: The first element is the actual entropy estimate
%%%  while the remaining 20 elements correspond to 20 distinct bootstrap
%%%  estimates (see Fig. 2).
%%%
%%%
%%%                Actual | Bootstrap
%%%              Estimate | Estimates
%%%                       |
%%%           Index:    1 | 2   3   4      Opt.btsp+1 = 21
%%%                   ----|------------- ... -----
%%%                   | x | x | x | x |      | x |
%%%                   ----|------------- ... -----
%%%                       |
%%%
%%% ### Building and calling custom bias correction functions
%%%
%%% In the direct method each value, P(x), of a probability distribution is
%%% estimated using the normalized count, C(x), of occurrence of x as 
%%% $P(x) = C(x)/N$.
%%%
%%% C contains all the information regarding the probability distribution
%%% P(x) and thus also all the parameters necessary for estimation of the
%%% bias of H(X). The same applies to the estimation of the other Y-
%%% unconditional probabilities, such as P_sh(x), and also to the Y-
%%% conditional probabilities: in the Y-conditional case bias
%%% correction is performed for each of the Ny probability distributions
%%% P(x|y) and then averaged across the values of Y.
%%%
%%% Users can define their custom bias correction methods as follows: let
%%% `custom_bias_corr_func` be the name of the user-defined correction routine
%%% (any valid MATLAB function name can be used, see Matlab documentation).
%%% This function must receive as input the array C of size *opts.n_binsX x
%%% 1* described above and return the (positive) value to be **added** to the
%%% plugin entropy estimate:
%%%       
%%%         bias = custom_bias_corr_func(C)
%%%
%%% To call the custom bias correction functon simply pass the name of the
%%% function as a string in the Opt.bias field. For the above example we
%%% have
%%%       
%%%       opts.bias = 'custom_bias_corr_func';
%%%
%%% ### Remarks
%%% - Field-names in the option structure are case-sensitive
%%% - Ouput options are case INsensitive
%%% - It is more efficient to call INFORMATION with several output options rather than calling the function repeatedly. For example:
%%% ```
%%% [X, Y] = information(X, Y, opts, {'I', 'Ish'});
%%% ```
%%% is faster than
%%% ```
%%% X = information(X, Y, opts, {'I'});
%%% Y = information(X, Y, opts, {'Ish'});
%%%```
%%%
