clc, clear;
% -------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox -  TE/cTE
% -------------------------------------------------------------------------
% -------------------------------------------------------------------------
%                                 TE
% -------------------------------------------------------------------------
% First we generate random data. The TE function TE() in the MINT Toolbox 
% works with input data in form of a cell. Suppose we have data from two neural populations X1, X2 and
% stimulus data S:
X1 = rand(2, 30, 100);
X2 = rand(2, 30, 100);

% MINT offers two TE outputs (you can also ask for both):
% - 'TE(A->B)' : Transfer entropy from A to B
% - 'TE(B->A)' : Transfer entropy from B to A 

% If we want to compute the transfer entropy from X1 to X2 and from X2 to X1, 
% we can use the following as input and output list:
inputs = {X1, X2};
outputList = {'TE(A->B)', 'TE(B->A)'};

% The MINT toolbox allows you to define optional parameters using a 
% structure with various fields. For TE computation, you can specify 
% the following options. If not provided, the function will use default 
% options automatically.

% As TE computes the transfer of entropy, we have to define the timepoints 
% in which we want to measure the transfer.
% Therefore we have to specify the timepoint which should be considered as
% 'present' and tau which defines the delays that should be used to compute
% the transfer of information. The parameters have to be defined in the
% variables opts.tpres and opts.tau. You can either just define one value,
% which will be then used for all Sources or you can define specific values
% for each source.
TE_opts.tau = {[3,5,15], [5, 7, 20]};      % Options:                                                      | (default: 1)
                                           % Integer or array with value(s) < tpres

TE_opts.tpres = {29, 30};                  % Integer that is > tau                                         | (default: last Timepoint of source A)

% The TE function offers two ways to compute transfer entropy based on the
% specified delays. Either it uses the data of that single timepoint or it
% uses timeseries delayed by the specified tau. Use opts.singleTimepoint to choose.
TE_opts.singleTimepoint = true;                    % Options:                                              | (default: false)
                                                   % true/false

% To bin your data, you must specify both the binning method and the 
% number of bins. You can set different methods or number of bins for each input; 
% if only one is provided, it will be applied to all inputs. 
% For more details type 'help binning' in your Command Window.
TE_opts.bin_method = {'eqpop', 'eqpop'};           % Options:                                              | (default: 'none')
                                                   % 'eqpop' (equal population), 
                                                   % 'eqspace' (equal spacing)
                                                   % 'threshold'
                                                   % 'none' (no binning)                                   
TE_opts.n_bins = {5, 5};                           % Specify an integer for the number of bins             | (default: 3)   


% To correct for sampling bias, you can specify methods for bias correction 
% in the MINT toolbox. If its not naive you will get the naive values as
% the  
% second output of the functions. 
% For more details type 'help correction' in your Command Window.
% You can define the bias correction options using the following parameters:
TE_opts.bias = 'shuffSub';                 % Options:                                                      | (default: 'naive')
                                            % 'naive' (no bias correction)                                                                   
                                            % 'qe' (Quadratic Extrapolation), 
                                            % 'le' (Linear Extrapolation), 
                                            % 'qe_shuffSub' (Quadratic with Shuffle Subtraction), 
                                            % 'le_shuffSub' (Linear with Shuffle Subtraction), 
                                            % 'shuffSub' (Shuffle Subtraction), 
                                            % 'pt' (Panzeri-Treves), 
                                            % 'bub' (BUB Correction), 
                                            % 'shuffCorr' (Shuffle Correction) 
% If you use 'qe' you have to specify the number of extrapolations (opts.xtrp) the
% correction function should to. (default: 10).
% If you use 'shuffSub' you have to specify the number of shufflings (opts.shuff) the
% correction function should to. (default: 20).
TE_opts.shuff = 20;   


% Several functions give you warning, e.g. if you did not specify an opts
% field that is needed and the function is using the default ist will
% inform you. If you don't want to get these warning you can supress them
% with the opts field supressWarnings (default: false).
TE_opts.supressWarnings = false;

% To test for significance, you may want to generate a null distribution.
% The functions in the MINT toolbox offer an option to do this through 
% the `opts.computeNulldist` parameter. To obtain the null distribution, 
% set this option to `true` and specify the number of samples using 
% `opts.n_samples` (default is 100).
%
% You can also indicate which variables should be shuffled to create the 
% null distribution or specify if the shuffling should be conditioned on 
% another variable using `opts.shuffling`. If you provide multiple 
% shuffling definitions, the function will output one null distribution 
% for each shuffling definition specified. (default: 'A' (just the first
% input variable is going to be shuffled). You can also specify the
% dimension that should be shuffled with opts.dim_shuffle. For 'Trials' it
% shuffles the last dimension, for 'Objects' the first and for 'Timepoints'
% the second.
% For more information, type 'help create_NullDistribution' in your Command 
% Window.
TE_opts.computeNulldist = true;
TE_opts.n_samples = 20;
TE_opts.shuffling = {'AB'};
TE_opts.dim_shuffle = {'Trials'};

% No we defined all options and we can call the TE function H() as
% follows: 
[TE_corrected, TE_naive, TE_nullDist] =   TE(inputs, outputList, TE_opts); 

% -------------------------------------------------------------------------
%                                 cTE
% -------------------------------------------------------------------------
% Conditioned transfer entropy between two time series conditioned 
% on a third time series. Therefore we need a third data, that the
% transferentropy is conditioned on:
X3 = rand(2, 30, 100);
inputs = {X1, X2, X3};
outputList = {'TE(A->B|C)', 'TE(B->A|C)'};

% Given the third dataset we also have to specify the opts for this data. 
% If we don't do that, the function would use the last specified bin_method,
% n_bins, tau and tpres for X3 (C). 
TE_opts.tau = {[3,5,15], [5, 7, 20], [5,8]};
TE_opts.tpres = {29, 30, 25};
TE_opts.bin_method = {'eqpop', 'eqpop', 'eqpop'}; 
TE_opts.n_bins = {5, 5, 6}; 

[cTE_corrected, cTE_naive, cTE_nullDist] = cTE(inputs, outputList, TE_opts); 









