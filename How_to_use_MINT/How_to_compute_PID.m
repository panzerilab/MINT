clc, clear;
% -------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox -  PID
% -------------------------------------------------------------------------
% The beginning of this Tutorial is a small simulation to generate data, on
% which we can show how to use the PID function of MINT and which options
% it offers.

% We simulate spiking activity of two neurons to four stimuli values as a
% sum of two poisson processes (one independent and one shared) for each
% neuron. 

% First we define the following parameters (you can change them to see how they influence the PID atoms):
B = 5;              % Baseline parameter (increasing B decreases the joint information) 
beta = 0.1;         % Dissimilarity of spike tuning (increasing beta increases redundancy)
gamma_flat = 20;    % strength of the shared process (increasing gamma increases the synergy)
alpha=10;           % Strength of stimulus tuning of the spike counts (increasing alpha increases the joint information)

% Based on the above specified parameter we now define the spike rates of
% the two neurons
rateX1 = B + alpha*beta*[1,1,0,0] + alpha*(1-beta)*[1,0,1,0];
rateX2 = B + alpha*beta*[1,0,1,0] + alpha*(1-beta)*[1,1,0,0];
rateShared = gamma_flat*[1,1,1,1];

num_stimuli = 4;
num_trials = 400;

% Now we simulate the spike trains for X1 and X2 as a sum of the individualspikes and the shared process and the stimulus S
X1 = [poissrnd(rateX1(1),1,num_trials/4) poissrnd(rateX1(2),1,num_trials/4)  poissrnd(rateX1(3),1,num_trials/4) poissrnd(rateX1(4),1,num_trials/4)];
X2 = [poissrnd(rateX2(1),1,num_trials/4) poissrnd(rateX2(2),1,num_trials/4)  poissrnd(rateX2(3),1,num_trials/4) poissrnd(rateX2(4),1,num_trials/4)];
Xshared = [poissrnd(rateShared(1),1,num_trials/4) poissrnd(rateShared(2),1,num_trials/4)  poissrnd(rateShared(3),1,num_trials/4) poissrnd(rateShared(4),1,num_trials/4)];
X1 = X1+Xshared;
X2 = X2+Xshared;
S = [ones(1,num_trials/4) 2*ones(1,num_trials/4) 3*ones(1,num_trials/4) 4*ones(1,num_trials/4)];

%% How to compute PID with MINT

% MINT offers the following PID outputs:
% - 'Syn'      : the synergistic component of the source neurons about the target
% - 'Red'      : the redundant component of the source neurons about the target
% - 'Unq1'     : the unique component of the first source neuron about the target
% - 'Unq2'     : the unique component of the second source neuron about the target
% - 'Unq'      : the sum of the unique components of the source neurons about the target
% - 'PID_atoms': All the PID atoms (for two sources: Syn, Red, Unq1, Unq2 [in that order])
% - 'Joint'    : The sum of all the PID atoms, equal to the joint information of the sources about the target 
% - 'Union'    : The sum of the atoms 'Red', 'Unq' 
% - 'q_dist'   : This output option is available only when the redundancy measure is set to 'I_BROJA'. 
%                The 'q_dist' provides the probability distribution derived from the Broja optimization 
%                process. This optimization seeks to maximize the joint information while maintaining the 
%                pairwise marginals at specified target values.
% - 'all'      : outputs all ouputs described above in the following order:
%                {'Syn', 'Red', 'Unq1', 'Unq2', 'Unq', 'Joint', 'Union'} (+ 'q_dist' for I_BROJA)

% If we want to compute the PID Atoms for X1 and X2 with S as target, we
% can use the following as input and output list:
inputs = {X1, X2, S};
outputList = {'PID_atoms'};

% The MINT toolbox allows you to define optional parameters using a 
% structure with various fields. For PID computation, you can specify 
% the following options. If not provided, the function will use default 
% options automatically.

% In the MINT toolbox you have different options for the redundancy measures that can be
% specified with the opts.redundancy_measure input option. 
PID_opts.redundancy_measure = 'I_min';      % Options:                                                  | (default: 'I_BROJA')
                                            % 'I_BROJA'
                                            % 'I_MMI'
                                            % 'I_min'                                              

% To bin your data, you must specify both the binning method and the 
% number of bins. You can set different methods or number of bins for each input; 
% If only one is provided, it will be applied to all inputs. 
% For more details type 'help binning' in your Command Window.
PID_opts.bin_method = {'eqpop', 'eqpop', 'none'};   % Options:                                          | (default: 'none')
                                                    % 'eqpop' (equal population), 
                                                    % 'eqspace' (equal spacing)
                                                    % 'threshold'
                                                    % 'none' (no binning)                                   
PID_opts.n_bins = {3, 3};                           % Specify an integer for the number of bins         | (default: 3)   

% To correct for limited-sampling bias, you can specify methods for bias correction 
% in the MINT toolbox. If its not naive you will get the naive values as
% the second output of the functions. 
% For more details type 'help correction' in your Command Window.
PID_opts.bias = 'naive';                            % Options:                                          | (default: 'naive')
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
PID_opts.pid_constrained = true;                    % Options.                                          | (default: true)
                                                    % true/false

% Several functions give you warning, e.g. if you did not specify an opts
% field that is needed and the function is using the default ist will
% inform you. If you don't want to get these warning you can supress them
% with the opts field supressWarnings (default: false).
PID_opts.supressWarnings = false;

% If you have real data, sometimes it includes NaN values. In the MINT
% toolbox you can specify how these should be handled with the field NaN
% handling. There are to options, either 'removeTrial' so these trials will
% be removed from all variables or 'error', so the compuation stops and it
% will give a warning that the data included NaN values.
PID_opts.NaN_handling = 'removeTrial';              % Options:                                          | (default: 'error')
                                                    % 'error' (Throws an error if NaN values are detected in any input)
                                                    % 'removeTrial' (Removes trials (across all variables) that contain NaN values)  

% To test for significance, you may want to generate a null distribution.
% The functions in the MINT toolbox offer an option to do this through 
% the `opts.computeNulldist` parameter. To obtain the null distribution, 
% set this option to `true` and specify the number of samples using 
% `opts.n_samples` (default is 100).
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
PID_opts.computeNulldist = false;  
PID_opts.parallel_sampling = true;   
PID_opts.n_samples = 20;
PID_opts.shuffling = {'AB', 'A'};
PID_opts.dim_shuffle = {'Trials'};

% No we defined all options and we can call the PID function PID() as follows: 
[PID_corrected, PID_naive, PID_nullDist] =   PID(inputs, outputList, PID_opts); 

