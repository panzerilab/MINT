clc, clear;
% -------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox -  PID
% -------------------------------------------------------------------------
% First we generate random data. The PID function PID() in the MINT Toolbox 
% works with input data in form of a cell. Suppose we have data from two neural populations X1, X2 and
% stimulus data S:
B = 5; beta = 0.1;gamma_flat = 20;alpha=10;
rate1 = B + alpha*beta*[1,1,0,0] + alpha*(1-beta)*[1,0,1,0];
rate2 = B + alpha*beta*[1,0,1,0] + alpha*(1-beta)*[1,1,0,0];
rateShared = gamma_flat*[1,1,1,1];
num_stimuli = length(rate1);
num_trials = 100;
X1=[poissrnd(rate1(1),1,num_trials) poissrnd(rate1(2),1,num_trials)  poissrnd(rate1(3),1,num_trials) poissrnd(rate1(4),1,num_trials)];
X2=[poissrnd(rate2(1),1,num_trials) poissrnd(rate2(2),1,num_trials)  poissrnd(rate2(3),1,num_trials) poissrnd(rate2(4),1,num_trials)];
Xshared=[poissrnd(rateShared(1),1,num_trials) poissrnd(rateShared(2),1,num_trials)  poissrnd(rateShared(3),1,num_trials) poissrnd(rateShared(4),1,num_trials)];
X1 = X1+Xshared;
X2 = X2+Xshared;
S=[ones(1,num_trials) 2*ones(1,num_trials) 3*ones(1,num_trials) 4*ones(1,num_trials)];

% -------------------------------------------------------------------------
%                                 PID
% -------------------------------------------------------------------------
% MINT offers various PID outputs:
% - 'Syn'      : the synergistic component of the source neurons about the target
% - 'Red'      : the redundant component of the source neurons about the target
% - 'Unq1'     : the unique component of the first source neuron about the target
% - 'Unq2'     : the unique component of the second source neuron about the target
% - 'Unq'      : the sum of the unique components of the source neurons about the target
% - 'PID_atoms': All the PID atoms (Syn, Red, Unq1, Unq2 [in that order])
% - 'Joint'    : The sum of all the PID atoms, equal to the joint information of the sources about the target 
% - 'Union'    : The sum of the atoms 'Red', 'Unq' 
% - 'q_dist'   : The probability distribution used to calculate the PID terms.

% - 'all'      :  outputs all ouputs described above in the following order: {'Syn', 'Red', 'Unq1', 'Unq2', 'Unq', 'Joint', 'Union', 'q_dist'}

% If we want to compute the PID Atoms for X1 and X2 with S as target, we
% can use the following as input and output list:
inputs = {X1, X2, S};
outputList = {'PID_atoms'};

% The MINT toolbox allows you to define optional parameters using a 
% structure with various fields. For PID computation, you can specify 
% the following options. If not provided, the function will use default 
% options automatically.

% There are different methods to decompose information. In the MINT
% toolbox you have different option for the redundancy measures that can be
% specified with the opts.redundancy_measure input option. 
PID_opts.redundancy_measure = 'I_BROJA';    % Options:                                                  | (default: 'I_BROJA')
                                            % 'I_BROJA'
                                            % 'I_MMI'
                                            % 'I_min'                                              

% To bin your data, you must specify both the binning method and the 
% number of bins. You can set different methods or number of bins for each input; 
% if only one is provided, it will be applied to all inputs. 
% For more details type 'help binning' in your Command Window.
PID_opts.bin_method = {'eqpop', 'eqpop', 'none'};   % Options:                                          | (default: 'none')
                                                    % 'eqpop' (equal population), 
                                                    % 'eqspace' (equal spacing)
                                                    % 'threshold'
                                                    % 'none' (no binning)                                   
PID_opts.n_bins = {3, 3};                           % Specify an integer for the number of bins         | (default: 3)   


% To correct for sampling bias, you can specify methods for bias correction 
% in the MINT toolbox. If its not naive you will get the naive values as
% the  
% second output of the functions. 
% For more details type 'help correction' in your Command Window.
% You can define the bias correction options using the following parameters:
PID_opts.bias = 'abc';                 % Options:                                                        | (default: 'naive')
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
PID_opts.xtrp = 10;   

% In the case of two sources and one target MINT offers the option for shuffSub and qe, to
% compute the corrected PID atoms constrained to the corrected information.
% This can be set with the option opts.pid_constrained:
PID_opts.pid_constrained = true;                % Options.                                              | (default: true)
                                            % true/false

% Several functions give you warning, e.g. if you did not specify an opts
% field that is needed and the function is using the default ist will
% inform you. If you don't want to get these warning you can supress them
% with the opts field supressWarnings (default: false).
PID_opts.supressWarnings = false;

PID_opts.NaN_handling = 'removeTrial';      % Options:                                                  | (default: 'error')
                                            % 'error' (Throws an error if NaN values are detected in any input)
                                            % 'removeTrial' (Removes trials (across all variables) that contain NaN values)  

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
PID_opts.computeNulldist = true;
PID_opts.n_samples = 20;
PID_opts.shuffling = {'AB', 'A'};
PID_opts.dim_shuffle = {'Trials'};

% No we defined all options and we can call the PID function PID() as follows: 
[PID_corrected, PID_naive, PID_nullDist] =   PID(inputs, outputList, PID_opts); 

