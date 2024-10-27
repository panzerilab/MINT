clc, clear;
% -------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox -  Mutual Information
% -------------------------------------------------------------------------
% First we generate random data. The Mutual Information function MI() in the MINT Toolbox 
% works with input data in form of a cell. Suppose we have neural data of
% X1 that encodesa stimulus S from timepoint 10 to 20:
rate1 = [12, 12, 2, 2];  
rate2 = [4, 4, 4, 4]; 
num_trials = 400;      
num_timepoints = 15;
stim_timepoints = 5:10;
X1 = zeros(1,num_timepoints, num_trials);
for tP = stim_timepoints
    X1(:, tP,:)=[poissrnd(rate1(1), 1, num_trials/4), ...
                 poissrnd(rate1(2), 1, num_trials/4), ...
                 poissrnd(rate1(3), 1, num_trials/4), ...
                 poissrnd(rate1(4), 1, num_trials/4)];

end
non_stim_timepoints = [1:4, 11:15];
for tP = non_stim_timepoints
    X1(:, tP,:) = [poissrnd(rate2(1), 1, num_trials/4), ...
                   poissrnd(rate2(2), 1, num_trials/4), ...
                   poissrnd(rate2(3), 1, num_trials/4), ...
                   poissrnd(rate2(4), 1, num_trials/4)];
end
S = [ones(1, num_trials/4), 2*ones(1, num_trials/4), 3*ones(1, num_trials/4), 4*ones(1, num_trials/4)];

% -------------------------------------------------------------------------
%                                 Mutual Information
% -------------------------------------------------------------------------
% MINT offers various Information measures:
% - 'I(A;B)'    : Mutual Information I(A;B)
% - 'Ilin(A;B)' : Linear MI Ilin(A;B)
% - 'coI(A;B)'  : Co-information coI(A;B)
% - 'Iss(A)'    : Sig. Sim Information Iss(A)
% - 'Ic(A;B)'   : Sum of Ici and Icd - Ic(A;B)
% - 'Ici(A;B)'  : Correlation Independent Information Ici(A;B)
% - 'Icd(A;B)'  : Correlation Dependent Information Icd(A;B)

% If we want to compute for example I(A;B), Ilin(A;B) and Ici(A;B) for X1 and X2; we specify
% the input and the outputs list as follows:
inputs = {X1, S};
outputList = {'I(A;B)', 'Ilin(A;B)', 'Ici(A;B)'};

% The MINT toolbox allows you to define optional parameters using a 
% structure with various fields. For MI computation, you can specify 
% the following options. If not provided, the function will use default 
% options automatically.

% To bin your data, you must specify both the binning method and the 
% number of bins. You can set different methods or number of bins for each input; 
% if only one is provided, it will be applied to all inputs. 
% For more details type 'help binning' in your Command Window.
MI_opts.bin_method = {'eqpop', 'none'};     % Options:                                      | (default: 'none')
                                            % 'eqpop' (equal population), 
                                            % 'eqspace' (equal spacing)
                                            % 'userEdges'
                                            % 'none' (no binning)                                   
MI_opts.n_bins = {4};                       % Specify an integer for the number of bins     | (default: 3)   


% To correct for sampling bias, you can specify methods for bias correction 
% in the MINT toolbox. If its not naive you will get the naive values as the 
% second output of the functions. 
% For more details type 'help correction' in your Command Window.
% You can define the bias correction options using the following parameters:
MI_opts.bias = 'shuffSub';             % Options:                                           | (default: 'naive')
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
MI_opts.shuff = 30;   

% Several functions give you warning, e.g. if you did not specify an opts
% field that is needed and the function is using the default ist will
% inform you. If you don't want to get these warning you can supress them
% with the opts field supressWarnings (default: false).
MI_opts.supressWarnings = true;

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
MI_opts.computeNulldist = true;
MI_opts.n_samples = 10;
MI_opts.shuffling = {'AB'};
MI_opts.dim_shuffle = {'Trials'};

% No we defined all options and we can call the MI function MI() as
% follows: 
[MI_corrected, MI_naive, MI_nullDist] =   MI(inputs, outputList, MI_opts); 

