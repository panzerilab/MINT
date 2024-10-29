clc, clear;
% -------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox -  Entropy
% -------------------------------------------------------------------------
% First we generate random data. The entropy function H() in the MINT Toolbox 
% works with input data in form of a cell. Suppose we have neural data from X1, X2, X3 and
% stimulus data S:
X1 = rand(2,10, 100);
X2 = rand(1,10, 100);
% -------------------------------------------------------------------------
%                                 Entropy
% -------------------------------------------------------------------------
% MINT offers various entropy measures:
% - 'H(A)'        : Measures the total uncertainty (entropy) of A.
% - 'H(A|B)'      : Entropy of A conditioned on B, representing uncertainty in A given B.
% - 'Hlin(A)'     : Linear approximation of the entropy of A.
% - 'Hind(A)'     : Independent entropy of A, assuming no correlations.
% - 'Hind(A|B)'   : Conditional independent entropy of A given B.
% - 'Chi(A)'      : Chi-squared-based entropy of A.
% - 'Hsh(A)'      : Entropy of A after shuffling.
% - 'Hsh(A|B)'    : Shuffled conditional entropy of A given B.

% If we want to compute for example H(A), H(A|B) and Hind(A|B) for X1 and X2 we specify
% the input and the outputs list as follows:
inputs = {X1, X2};
%outputList = {'H(A)', 'H(A|B)', 'Hind(A|B)'};
outputList = {'H(A)', 'H(A|B)'};

% The MINT toolbox allows you to define optional parameters using a 
% structure with various fields. For entropy computation, you can specify 
% the following options. If not provided, the function will use default 
% options automatically.

% To bin your data, you must specify both the binning method and the 
% number of bins. You can set different methods or number of bins for each input; 
% if only one is provided, it will be applied to all inputs. 
% For more details type 'help binning' in your Command Window.
entropy_opts.bin_method = {'eqpop', 'eqspace'};     % Options:                                      | (default: 'none')
                                                    % 'eqpop' (equal population), 
                                                    % 'eqspace' (equal spacing)
                                                    % 'userEdges
                                                    % 'none' (no binning)                                   
entropy_opts.n_bins = {4, 3};               % Specify an integer for the number of bins             | (default: 3)   


% To define how NaN values in your data should be handled, you can specify 
% opts.NaN_handling. 
entropy_opts.NaN_handling = 'removeTrial';  % Options:                                              | (default: 'error')
                                            % 'error' (Throws an error if NaN values are detected in any input)
                                            % 'removeTrial' (Removes trials (across all variables) that contain NaN values)                                             
                                            % 'setToZero' (Replaces NaN values with zeros in all input data)                                                           


% To correct for sampling bias, you can specify methods for bias correction 
% in the MINT toolbox. If its not naive you will get the naive values as the 
% second output of the functions. 
% For more details type 'help correction' in your Command Window.
% You can define the bias correction options using the following parameters:
entropy_opts.bias = 'bub';             % Options:                                              | (default: 'naive')
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
entropy_opts.shuff = 20;   

% Several functions give you warning, e.g. if you did not specify an opts
% field that is needed and the function is using the default ist will
% inform you. If you don't want to get these warning you can supress them
% with the opts field supressWarnings (default: false).
entropy_opts.supressWarnings = false;

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
entropy_opts.computeNulldist = true;
entropy_opts.n_samples = 100;
entropy_opts.parallel_sampling = false;
entropy_opts.shuffling = {'AB','A'};
entropy_opts.dim_shuffle = {'Trials'};

% No we defined all options and we can call the entropy function H() as
% follows: 
[entropy_corrected, entropy_naive, entropy_nullDist] =   H(inputs, outputList, entropy_opts); 


