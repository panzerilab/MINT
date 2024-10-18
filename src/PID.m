function [PID_values, PID_naive, PID_shuff_all] = PID(inputs, varargin)
%%% PID - Calculate Partial Information Decomposition (PID) and related information-theoretic quantities
%%%
%%% This function calculates the atoms partial information decomposition 
%%% (PID) and other related measures based on the provided inputs, outputs,
%%% and optional parameters.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input time series data. Each cell represents a time series, where:
%%%             - inputs{1}: First time series (A) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%             - inputs{2}: Second time series (B) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%   - outputs: A cell array of strings specifying which entropies to compute.
%%%              Possible outputs include:
%%%               - 'Syn'    : the synergistic component of the source neurons about the target(2 sources case)
%%%               - 'Red'    : the redundant component of the source neurons about the target(2 sources case)
%%%               - 'Unq1'   : the unique component of the first source neuron about the target(2 sources case)
%%%               - 'Unq2'   : the unique component of the second source neuron about the target(2 sources case)
%%%               - 'Unq'    : the sum of the unique components of the source neurons about the target(2 sources case)
%%%               - 'PID_atoms'  : All the PID atoms (n sources case) (default)
%%%               - 'Joint'  : The sum of all the PID atoms, equal to the joint information of the sources about the target (n sources case)
%%%               - 'Union'  : The sum of the atoms 'Red', 'Unq' (2 sources case)
%%%               - 'q_dist' : The probability distributionused to calculate the PID terms.
%%%
%%%   - varargin: Optional arguments, passed as a structure. Fields may include:
%%%              - redundancy_measure: name of the measure of the redundancy between sources 
%%%                                Possible values include:
%%%                                'I_BROJA'                    :(default)
%%%                                only available for two sources, which
%%%                                defines redundancy as the result of a
%%%                                constrained optimization problem (Bertschinger et al., 2014; Makkeh et al., 2018)
%%%                                'I_MMI'                      : minimal
%%%                                 mutual information (MMI) defines the 
%%%                                 redundancy as the smallest single
%%%                                 information between a source and the
%%%                                 target (Barret, 2015)
%%%                                'I_min'                      : redundancy measure proposed by (Williams and Beer, 2010)
%%%                                
%%%                                Users can also define their own custom bias correction method
%%%                                (see help for correction.m)
%%%
%%%              - bias:           Specifies the bias correction method to be used.
%%%                                Possible values include:
%%%                                'naive'                      :(default) - No correction applied.
%%%                                'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%%%                                'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%%%                                'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).

%%%                                Users can also define their own custom bias correction method
%%%                                (see help for correction.m)
%%%
%%%              - bin_method:     Cell array specifying the binning method.
%%%                                It can have one or two entries:
%%%                                If one entry is provided, it will be applied to both A and B.
%%%                                Possible values include:
%%%                                'eqpop'     : Equal population binning.
%%%                                'eqspace'   : Equal space binning.
%%%                                'threshold' : Binning based on a specified threshold.
%%%                                Users can also define their own custom binning method
%%%                                (see help for binning.m).
%%%                                Default is {'none'}.
%%%
%%%              - n_bins:         Specifies the number of bins to use for binning.
%%%                                It can be a single integer or a cell array with one or two entries.
%%%                                If one entry is provided, it will be used for both A and B.
%%%                                This integer defines how the continuous values will be
%%%                                discretized into bins for analysis.
%%%                                Default number of bins is {3}.
%%%
%%%              - suppressWarnings: Boolean (true/false) to suppress warning messages.
%%%                                  Default is false, meaning warnings will be shown.
%%% Outputs:
%%%   - PID_values: A cell array containing the computed MI values as specified in the outputs argument.
%%%   - PID_naive: A cell array containing the naive MI estimates.
%%%   - PID_shuff_all: A value indicating the all results of the shuffling procedure (0 if not performed).
%%%
%%% EXAMPLE
%%% Suppose we have two groups of neurons X1 and X2
%%% and a group of neurons Y.
%%% (Structure of X1, X2, Y is nNeurons x nTrials)
%%%
%%% We can structure our inputs as follows:
%%% Thus, the total input for A and B would be:
%%% A = {X1; X2}; 
%%% B = Y;  % 
%%%
%%% To compute the synergy and redundancy between the sources (A) about the target (B), the function can be called as:
%%% PID_values = PID({A, B}, {'Syn', 'Red'}, opts);
%%%
%%% Here, 'opts' represents additional options you may want to include, such as
%%% specifying the bias correction method, number of bins (n_bins), and other parameters as needed.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check Inputs, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            outputs = varargin{1};
        end
    else
        [outputs, opts] = check_inputs('PID',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('PID',inputs,varargin{:});
end

possibleOutputs = {'', 'Syn',  'Red', 'Unq1', 'Unq2', ...
    'Unq', 'PID_atoms', 'Joint', 'Union', 'q_dist'};
if ismember('', outputs)
    outputs = {'Syn', 'Red', 'Unq1', 'Unq2'};
elseif ismember('PID_atoms', outputs)
    outputs(ismember(outputs, 'PID_atoms')) = [];  
    outputs = [outputs, {'Syn', 'Red', 'Unq1', 'Unq2'}];
elseif ismember('all', outputs)
    outputs = {'Syn', 'Red', 'Unq1', 'Unq2', 'Unq', 'Joint', 'Union', 'q_dist'};
end
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('PID:invalidOutput', msg);
end

nSources = length(inputs)-1;
if nSources < 2
    msg = 'At least two sources are required.';
    error('PID:NotEnoughSources', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Binning, reduce dimensions if necessary                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~opts.isBinned == false
    inputs_b = binning(inputs, opts);   
    opts.isBinned = true;
else
    inputs_b = inputs;
end

inputs_1d = inputs_b;
nTrials_comp = size(inputs_1d{1},2);
for var = 1:nSources+1
   sizeVar = size(inputs_1d{var});
   nTrials = sizeVar(end);
   if nTrials ~= nTrials_comp
       msg = 'Inconsistent input size. Number of Trials must be consistent.';
       error('PID:Invalid Input', msg);
   end 
   if sizeVar(1) > 1
       inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
   end 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 3.A: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr = opts.bias;
corefunc = @PID;
if ~strcmp(corr, 'naive')
    [PID_values, PID_naive, PID_shuff_all] = correction(inputs_1d, outputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 3.B: Compute required PID Atoms                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch opts.redundancy_measure
    case 'I_BROJA'
        if nSources == 2
            opts.function = @pidBROJA;
        else
            warning('BROJA can be used for two sources only, switching to I_min')
            opts.redundancy_measure = 'I_min';
            opts.function = @pidimin;
        end
    case 'I_MMI'
        opts.function = @pidimmi;
    case 'I_min'
        opts.function = @pidimin;
end
opts.multidim = true;
sourceInput = cat(1, inputs_1d{1:end-1});
[p_distr] = prob_estimator({inputs_1d{end}, sourceInput}, {'P(A,B,C)'}, opts);
if strcmp(opts.redundancy_measure,'I_BROJA')
    [PID_terms, table_prob] = opts.function(p_distr);
else
    PID_terms = opts.function(p_distr);
end

if opts.pid_constrained
    I1 = cell2mat(MI({inputs_1d{1}, inputs_1d{end}}, {'I(A;B)'}, opts));
    I2  = cell2mat(MI({inputs_1d{2}, inputs_1d{end}}, {'I(A;B)'}, opts));
    I12 = cell2mat(MI({cat(1, inputs_1d{1}, inputs_1d{2}), inputs_1d{end}}, {'I(A;B)'}, opts));
    if ~isfield(opts, 'chosen_atom')
        opts.chosen_atom = 'red';
    end
    if strcmp(opts.chosen_atom, 'red')
        red = PID_terms(1);
        PID_terms = [red, I1-red, I2-red, I12-I1-I2+red];
    elseif strcmp(opts.chosen_atom, 'un1')
        un1 = PID_terms(2);
        PID_terms = [I1-un1, un1, I2-I1+un1, I12-I2-un1];
    elseif strcmp(opts.chosen_atom, 'un2')
        un2 = PID_terms(3);
        PID_terms = [I2-un2, I1-I2+un2, un2, I12-I1-un2];
    elseif strcmp(opts.chosen_atom, 'syn')
        syn = PID_terms(4);
        PID_terms = [I1+I2-I12+syn, I12-I2-syn, I12-I1-syn , syn];
    end
end


for i = 1:length(outputs) 
    output = outputs{i};
    switch output
        case 'Syn'
            PID_values{i} = PID_terms(4);
        case 'Red'
            PID_values{i} = PID_terms(1);
        case 'Unq1'
           PID_values{i} = PID_terms(2);
        case 'Unq2'
            PID_values{i} = PID_terms(3);
        case 'Unq'
            PID_values{i} = PID_terms(2)+PID_terms(3);
        case 'Joint'
            PID_values{i} = sum(PID_terms);
        case 'Union'
            PID_values{i} = sum(PID_terms)-PID_terms(4);
        case 'q_dist'
            if strcmp(opts.bias,'naive') && strcmp(opts.redundancy_measure,'I_BROJA') && istable(table_prob)
                q_distr = zeros(size(p_distr));
                for irow=1:height(table_prob)
                    x = table_prob{irow,'X'};
                    y = table_prob{irow,'Y'};
                    z = table_prob{irow,'Z'};
                    q_distr(x, y, z) = table_prob{irow, 'qpos'};
                end
                PID_values{i} = q_distr;
            else 
                PID_values{i} = NaN;
                if ~opts.supressWarnings
                    warning("q_dist can only be computed when the redundancy measure is set to 'I_Broja'.")
                end 
            end 
    end
end
PID_naive = PID_values;
end 
