clc, clear;
% -------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox -  II
% -------------------------------------------------------------------------
% First we generate random data. The II function II() in the MINT Toolbox 
% works with input data in form of a cell. Suppose we have data from a neural population X and
% stimulus data S as well as choice data C:
rate1 = [12, 12, 2, 2];  
rate2 = [4, 4, 4, 4]; 
num_trials = 400;      
num_timepoints = 15;
stim_timepoints = 5:10;
num_neurons = 2;
X = zeros(num_neurons,num_timepoints, num_trials);
for neuron = 1:num_neurons
    for tP = stim_timepoints
        X(neuron, tP,:)=[poissrnd(rate1(1), 1, num_trials/4), ...
                     poissrnd(rate1(2), 1, num_trials/4), ...
                     poissrnd(rate1(3), 1, num_trials/4), ...
                     poissrnd(rate1(4), 1, num_trials/4)];
    
    end
end 
non_stim_timepoints = [1:4, 11:15];
for neuron = 1:num_neurons
    for tP = non_stim_timepoints
        X(neuron, tP,:) = [poissrnd(rate2(1), 1, num_trials/4), ...
                       poissrnd(rate2(2), 1, num_trials/4), ...
                       poissrnd(rate2(3), 1, num_trials/4), ...
                       poissrnd(rate2(4), 1, num_trials/4)];
    end
end 
S = [ones(1, num_trials/4), 2*ones(1, num_trials/4), 3*ones(1, num_trials/4), 4*ones(1, num_trials/4)];
C = [1*ones(1, num_trials/4), 1*ones(1, num_trials/4), 2*ones(1, num_trials/4), 2*ones(1, num_trials/4)];

% -------------------------------------------------------------------------
%                                 II
% -------------------------------------------------------------------------
% MINT offers two II outputs (you can also ask for both):
% - 'II(A->B;S)' : Transfer entropy from A to B given S.
% - 'II(B->A;S)' : Transfer entropy from B to A given S.

% If we want to compute the feature information transfer from X1 to X2
% about S and from X2 to X1 about S, we can use the following as input and 
% output list:
inputs = {S, X, C};
outputList = {'II(A,B,C)'}; %You can also call the function without specifying the output here

% The MINT toolbox allows you to define optional parameters using a 
% structure with various fields. For II computation, you can specify 
% the following options. If not provided, the function will use default 
% options automatically.

% There are different methods to decompose information. In the MINT
% toolbox you have different option for the redundancy measures that can be
% specified with the opts.redundancy_measure input option. 
II_opts.redundancy_measure = 'I_BROJA';     % Options:                                               | (default: 'I_min')
                                            % 'I_MMI'
                                            % 'I_min' 
                                            % 'I_BROJA'

% To bin your data, you must specify both the binning method and the 
% number of bins. You can set different methods or number of bins for each input; 
% if only one is provided, it will be applied to all inputs. 
% For more details type 'help binning' in your Command Window.
II_opts.bin_method = {'none', 'eqpop', 'none'};   % Options:                                         | (default: 'none')
                                                    % 'eqpop' (equal population), 
                                                    % 'eqspace' (equal spacing)
                                                    % 'threshold'
                                                    % 'none' (no binning)                                   
II_opts.n_bins = {3};                               % Specify an integer for the number of bins      | (default: 3)   


% To correct for sampling bias, you can specify methods for bias correction 
% in the MINT toolbox. If its not naive you will get the naive values as
% the  
% second output of the functions. 
% For more details type 'help correction' in your Command Window.
% You can define the bias correction options using the following parameters:
II_opts.bias = 'qe';                 % Options:                                                      | (default: 'naive')
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
II_opts.xtrp = 10;   

% In the case of two sources and one target MINT offers the option for shuffSub and qe, to
% compute the corrected II atoms constrained to the corrected information.
% This can be set with the option opts.II_constrained:
II_opts.pid_constrained = true;            % Options.                                                 | (default: true)
                                            % true/false

% Several functions give you warning, e.g. if you did not specify an opts
% field that is needed and the function is using the default ist will
% inform you. If you don't want to get these warning you can supress them
% with the opts field supressWarnings (default: false).
II_opts.supressWarnings = false;

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
% For more informtion, type 'help create_nullDist' in your Command Window.
II_opts.computeNulldist = true;
II_opts.parallel_sampling = true;
II_opts.n_samples = 20;
II_opts.shuffling = {'ABC'};
II_opts.dim_shuffle = {'Trials'};

% No we defined all options and we can call the II function H() as
% follows: 
[II_corrected, II_naive, II_nullDist] = II(inputs, outputList, II_opts); 

% -------------------------------------------------------------------------
%                         tools: iiBoundaries
% -------------------------------------------------------------------------
opts.svm_options = struct();
opts.ii_options = struct();

% Make the stimulus binary, as iiBoundaries is a tool exclusive for binary stimulus and choice
S(S == 1 | S == 2) = 1;              
S(S == 3 | S == 4) = 2;           
inputs = {S, squeeze(X(:,10,:)), C};
outputs = iiBoundaries(inputs, {'all'}, opts);

% The outputs can include:
% - `II`: Intersection information between S R and C.
% - `theta`: Angle (in degrees) between the decision boundaries of S and C.
% - `labelsStim`: Predicted labels for Stimulus S by the SVM based on neural responses.
% - `labelsChoice`: Predicted labels for Choice C by the SVM based on neural responses
II_values    = outputs{1};                    
theta        = outputs{2};                     
labelsStim   = outputs{3}; 
labelsChoice = outputs{3}; 
                  








