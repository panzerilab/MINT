%%% *function [R_binned] = binr(R, nb, method, par)*
%%%
%%% ### Description
%%% Discretize the input matrix *R* .
%%%
%%% ### Inputs:
%%% - *R*: *nDims X nTrials* input matrix to be discretized
%%% - *nb*: number of bins array (can be a scalar or a vector of `size(nb) = size(R,1)`). If a scalar is specified then all dimensions are binned using the same number of bins. Alternatively, each dimension is binned using its specific number of bins. 
%%% - *method*: binning method (can be a scalar or a cell array of `size(method) = size(R,1)`). If a scalar is specified then all dimensions are binned using the same method. Alternatively, each dimension is binned using its specific method.
%%% - *par*: optional parameters
%%%
%%% ### Outputs:
%%% - *R_binned*: *nDims X nTrials* binned matrix.
%%%
%%% ### Further notes
%%% #### Binning methods and parameters
%%% The table below describe the binning options which are built-in in the toolbox. Alternatively, users can define and plug-in their custom binning methods (see [Building and calling custom binning functions](Building and calling custom binning functions)"Building and calling custom binning functions" section).
%%%
%%% Some of the built-in binning options allow a parameter to be specified.
%%% Parameters can be passed in a similar fashion to the options:
%%% - if a single parameter needs to be specified for all values of *R*, it can be passed directly or in the form of a one dimensional cell-array;
%%% - if different parameters need to be specified for the different binning options or values of *R*, then a cell array of length *nDims* must be passed to the function where each field specifies the binning option for the i-th value. To skip the parameter for a value, just pass an empty array [] for the corresponding index.
%%%
%%%   | Option        | Description                                           |
%%%   |---------------|-------------------------------------------------------|
%%%   | `'eqpop'`     | Equipopulated binning: the width of the bins is selected so that each bin contains (approximately) the same number of values. Using `'eqpop'` option for *R* values which are not continuous in nature (i.e., which contain several repeated values) may result in the bins being poorly equipopulated. |
%%%   | `'eqspace' `  | Equispaced binning: the range [a, b], provided as a parameter by the user in the form of a 2-element array, is divided into bins of equal width. If no interval is specified, then the range [m, M] (m and M being the max and min of the values to be binned, respectively) is used. A check is performed on whether the specified interval includes all values to be binned. |
%%%   | `'ceqspace'`  | Centered equispaced binning: the range [C-D,C+D], C being the mean of the values to be binned and D = max(C-m, M-C) (m and M being defined as above), is divided into intervals of equal width.
%%%   | `'gseqspace'` | Gaussian equispaced binning: the range [C-N*STD, C+N*STD] (STD being the standard deviation of the values to be binned, C being defined as above and N being a parameter which is passed to the function) is divided into intervals of equal width. If no N is specified, N=2 is used. If any of the values falls outside the selected range the first and last bin are stretched in order to accomodate outliers falling below or above the range limits, respectively.                           |
%%%
%%% #### Building and calling custom binning functions
%%% Users can define their custom binning methods by means of a binning routine of the form:
%%% ```      
%%% edges = <func_name>(x, nb)
%%% ```
%%% where:
%%% - `<func_name>`: any valid function name (see Matlab documentation for instructions on how to build valid function names)
%%% - `x`: a column array contatining the values which need to be quantized.
%%% - `nb`: the number of bin that need to be used for discretization.
%%% - `edges`: an `nb+1`-long array of strictly monotonically increasing values corresponding to the edges of the quantization bins
%%% To call the custom plug-in functons simply pass `'func_name'` as a string as the *method* parameter.
%%%
