clc, clear;
% -------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox -  FIT
% -------------------------------------------------------------------------
% First we generate random data. The FIT function FIT() in the MINT Toolbox 
% works with input data in form of a cell. Suppose we have data from two neural populations X1, X2 and
% stimulus data S:
rate1 = [15, 15, 5, 5]; 
rate2 = [4, 4, 4, 4]; 
num_trials = 200;      
num_timepoints = 15;
stim_timepoints = 1:10;
X1 = zeros(1,num_timepoints, num_trials);
X2 = zeros(1,num_timepoints, num_trials);
for tP = stim_timepoints
    X1(:, tP,:)=[poissrnd(rate1(1), 1, num_trials/4), ...
                 poissrnd(rate1(2), 1, num_trials/4), ...
                 poissrnd(rate1(3), 1, num_trials/4), ...
                 poissrnd(rate1(4), 1, num_trials/4)];

end
non_stim_timepoints = 11:15;
for tP = non_stim_timepoints
    X1(:, tP,:) = [poissrnd(rate2(1), 1, num_trials/4), ...
                   poissrnd(rate2(2), 1, num_trials/4), ...
                   poissrnd(rate2(3), 1, num_trials/4), ...
                   poissrnd(rate2(4), 1, num_trials/4)];
end
delay = 5;
for tP = (delay+1):num_timepoints
    X2(:, tP,:) = X1(:, tP-delay, :) + poissrnd(3,1,1,num_trials);
end
for tP = 1:delay
    X2(:, tP,:) = [poissrnd(rate2(1), 1, num_trials/4), ...
                   poissrnd(rate2(2), 1, num_trials/4), ...
                   poissrnd(rate2(3), 1, num_trials/4), ...
                   poissrnd(rate2(4), 1, num_trials/4)];
end


S = [ones(1, num_trials/4), 2*ones(1, num_trials/4), 3*ones(1, num_trials/4), 4*ones(1, num_trials/4)];
% -------------------------------------------------------------------------
%                                 FIT
% -------------------------------------------------------------------------
% MINT offers two FIT outputs (you can also ask for both):
% - 'FIT(A->B;S)' : Transfer entropy from A to B given S.
% - 'FIT(B->A;S)' : Transfer entropy from B to A given S.

% If we want to compute the feature information transfer from X1 to X2
% about S and from X2 to X1 about S, we can use the following as input and 
% output list:
inputs = {X1, X2, S};
outputList = {'FIT(A->B;S)', 'FIT(B->A;S)'};

% The MINT toolbox allows you to define optional parameters using a 
% structure with various fields. For FIT computation, you can specify 
% the following options. If not provided, the function will use default 
% options automatically.

% As FIT computes the transfer of information, we have to define the
% timepoints in which we want to measure the information transfer.
% Therefore we have to specify the timepoint which should be considered as
% 'present' and tau which defines the delays that should be used to compute
% the transfer of information. The parameters have to be defined in the
% variables opts.tpres and opts.tau. You can either just define one value,
% which will be then used for all Sources or you can define specific values
% for each source.
FIT_opts.tau = {5};                         % Options:                                                      | (default: 1)
                                            % Integer that is < tpres

FIT_opts.tpres = {15};                      % Integer that is > tau                                         | (default: last Timepoint of source A)

% There are different methods to decompose information. In the MINT
% toolbox you have different option for the redundancy measures that can be
% specified with the opts.redundancy_measure input option. 
FIT_opts.redundancy_measure = 'I_min';      % Options:                                                      | (default: 'I_min')
                                            % 'I_MMI'
                                            % 'I_min'                                              

% To bin your data, you must specify both the binning method and the 
% number of bins. You can set different methods or number of bins for each input; 
% if only one is provided, it will be applied to all inputs. 
% For more details type 'help binning' in your Command Window.
FIT_opts.bin_method = {'eqpop', 'eqpop', 'none'};   % Options:                                              | (default: 'none')
                                                    % 'eqpop' (equal population), 
                                                    % 'eqspace' (equal spacing)
                                                    % 'threshold'
                                                    % 'none' (no binning)                                   
FIT_opts.n_bins = {5, 5};                           % Specify an integer for the number of bins             | (default: 3)   


% To correct for sampling bias, you can specify methods for bias correction 
% in the MINT toolbox. If its not plugin you will get the plugin values as
% the  
% second output of the functions. 
% For more details type 'help correction' in your Command Window.
% You can define the bias correction options using the following parameters:
FIT_opts.bias = 'shuffSub';                 % Options:                                                      | (default: 'plugin')
                                            % 'plugin' (no bias correction)                                                                   
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
FIT_opts.xtrp = 10;   

% In the case of two sources and one target MINT offers the option for shuffSub and qe, to
% compute the corrected FIT atoms constrained to the corrected information.
% This can be set with the option opts.FIT_constrained:
FIT_opts.pid_constrained = true;            % Options.                                                       | (default: true)
                                            % true/false

% Several functions give you warning, e.g. if you did not specify an opts
% field that is needed and the function is using the default ist will
% inform you. If you don't want to get these warning you can supress them
% with the opts field supressWarnings (default: false).
FIT_opts.supressWarnings = true;

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
FIT_opts.computeNulldist = true;
FIT_opts.n_samples = 5;
FIT_opts.shuffling = {'AB_C'};
FIT_opts.dim_shuffle = {'Trials'};
FIT_opts.parallel_sampling = true;

% No we defined all options and we can call the FIT function H() as
% follows: 
%[FIT_corrected, FIT_plugin, FIT_nullDist] = FIT(inputs, outputList, FIT_opts); 


% -------------------------------------------------------------------------
%                                 cFIT
% -------------------------------------------------------------------------
% Conditioned transfer entropy between two time series conditioned 
% on a third time series. Therefore we need a third data, that the
% transferentropy is conditioned on:
X3 = zeros(1,num_timepoints, num_trials);
for tP = 1:5
    X3(:, tP,:)=[poissrnd(rate1(1), 1, num_trials/4), ...
                 poissrnd(rate1(2), 1, num_trials/4), ...
                 poissrnd(rate1(3), 1, num_trials/4), ...
                 poissrnd(rate1(4), 1, num_trials/4)];

end
for tP = 6:15
    X3(:, tP,:) = [poissrnd(rate2(1), 1, num_trials/4), ...
                   poissrnd(rate2(2), 1, num_trials/4), ...
                   poissrnd(rate2(3), 1, num_trials/4), ...
                   poissrnd(rate2(4), 1, num_trials/4)];
end
inputs = {X1, X2, X3, S};
outputList = {'cFIT(A->B;S|C)'};

% Given the third dataset we also have to specify the opts for this data. 
% If we don't do that, the function would use the last specified bin_method,
% n_bins, tau and tpres for X3 (C). 
FIT_opts.tau = {5};  
FIT_opts.tpres = {15};
FIT_opts.bin_method = {'eqpop', 'eqpop', 'eqpop','none'}; 
FIT_opts.n_bins = {5, 5, 5}; 
FIT_opts.shuffling = {'AB_D'};
FIT_opts.n_trials = num_trials*4;
%[cFIT_corrected, cFIT_plugin, cFIT_nullDist] = cFIT(inputs, outputList, FIT_opts); 
[cFIT_corrected, cFIT_plugin, cFIT_nullDist] = cFIT_old(S, X1,X2,X3, FIT_opts); 
