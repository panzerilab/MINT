function [PID_values, PID_naive, PID_nullDist] = PID(inputs, varargin)
% PID - Calculate Partial Information Decomposition (PID) and related information-theoretic quantities
%
% This function calculates the atoms partial information decomposition
% (PID) and other related measures based on the provided inputs, reqOutputs,
% and optional parameters.
%
% Inputs:
%   - inputs: A cell array containing the input data with N sources and one target:
%             - inputs{1}:   First source (A) with dimensions
%                            nDims X (nTimepoints X) nTrials
%             - inputs{2}:   Second source (B) with dimensions
%                            nDims X (nTimepoints X) nTrials
%                                    ...
%             - inputs{n+1}: Target data with dimensions
%                            nDims X (nTimepoints X) nTrials
%             -> In cases where the input is provided as a time series, the outputs 
%                will be computed for each time point, resulting in outputs that are 
%                also represented as time series.
%
%   - reqOutputs: A cell array of strings specifying which measures to compute:
%               - 'Syn'        : the synergistic component of the source neurons about the target
%               - 'Red'        : the redundant component of the source neurons about the target
%               - 'Unq1'       : the unique component of the first source neuron about the target
%               - 'Unq2'       : the unique component of the second source neuron about the target
%               - 'Unq'        : the sum of the unique components of the source neurons about the target
%               - 'PID_atoms'  : All the PID atoms 
%               - 'Joint'      : The sum of all the PID atoms, equal to the joint information of the sources about the target 
%               - 'Union'      : The sum of the atoms 'Red', 'Unq'
%               - 'q_dist'     : This output option is available only when the redundancy measure is set to 'I_BROJA'. 
%                                The 'q_dist' provides the probability distribution derived from the Broja optimization 
%                                process. This optimization seeks to maximize the joint information while maintaining the 
%                                pairwise marginals at specified target values.

%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - redundancy_measure: name of the measure of the redundancy between sources 
%                                    'I_BROJA' : only available for two sources, which defines redundancy as the result of a
%                                                constrained optimization problem (Bertschinger et al., 2014; Makkeh et al., 2018)
%                                    'I_MMI'   : minimal mutual information (MMI) defines the redundancy as the smallest 
%                                                singleinformation between a source and the target (Barret, 2015)
%                                    'I_min'   : redundancy measure proposed by (Williams and Beer, 2010)                              
%     
%              - bias:               Specifies the bias correction method to be used.
%                                    'naive'                      :(default) - No correction applied.
%                                    'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%                                    'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%                                    'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%                                    Users can also define their own custom bias correction method
%                                    (type 'help correction' for more information)
%     
%              - bin_method:         Cell array specifying the binning method to be applied.
%                                    'none'      : (default) - No binning applied.
%                                    'eqpop'     : Equal population binning.
%                                    'eqspace'   : Equal space binning.
%                                    'userEdges' : Binning based on a specified edged.
%                                    Users can also define their own custom binning method
%                                    If one entry is provided, it will be applied to both A and B.
%                                    (type 'help binning' for more information).
%     
%              - n_bins:             Specifies the number of bins to use for binning.
%                                    It can be a single integer or a cell array with one or two entries.
%                                    Default number of bins is {3}.
%
%              - computeNulldist:    If set to true, generates a null distribution
%                                    based on the specified inputs and core function.
%                                    When this option is enabled, the following can be specified:
%                                     - `n_samples`: The number of null samples to generate (default: 100).
%                                     - 'shuffling': Additional shuffling options to determine the variables to be 
%                                        shuffled during the computation of the null distribution (default: {'A'}).
%                                        (type 'help hShuffle' for more information).
%   
%              - suppressWarnings:    Boolean (true/false) to suppress warning messages.
%                                     Default is false, meaning warnings will be shown
%
%              - NaN_handling:     Specifies how NaN values should be handled in the data.
%                                  Options include:
%                                  'removeTrial' : Removes trials containing NaN in any variable 
%                                                  from all input data.
%                                  'error'       : (default) Throws an error if NaN values are detected.
%
%              - pid_constrained:   Boolean (true/false) to specify whether the partial information decomposition (PID) 
%                                   should be computed with constraints on specific atomic terms. When enabled, it uses 
%                                   the chosen atomic term specified in `chosen_atom` to compute the other atoms.
%                                   Default istrue, meaning PID is calculated with these constraints.
%                                   When this option is enabled, the following can be specified:
%                                    - `chosen_atom`: Determines the atomic term used in constrained PID calculations. 
%                                       Options include:
%                                       'Syn' : Synaptic or synergistic information (default).
%                                       'Red' : Redundant information.
%                                       'Unq1': Unique information from source 1.
%                                       'Unq2': Unique information from source 2.
%
%
% Outputs:
%   - PID_values: A cell array containing the computed MI values as specified in the reqOutputs argument.
%   - PID_naive: A cell array containing the naive MI estimates.
%   - PID_nullDist: Results of the null distribution computation (0 if not performed).
%
% EXAMPLE
% Suppose we have two groups of neurons X1 and X2 and a Stimulus S.
% To compute the synergy and redundancy between the sources X1 and X2 about 
% the target S, the function can be called as:
% PID_values = PID({X1, X2, S}, {'Syn', 'Red}, opts);


% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check Inputs, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 1
    msg = 'Please input your data.';
    error('PID:notEnoughInput', msg);
end

if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            reqOutputs = varargin{1};
        end
    else
        [inputs, reqOutputs, opts] = check_inputs('PID',inputs,varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('PID',inputs,varargin{:});
end

possibleOutputs = {'', 'Syn',  'Red', 'Unq1', 'Unq2', ...
    'Unq', 'PID_atoms', 'Joint', 'Union', 'q_dist'};
if ismember('', reqOutputs)
    reqOutputs = {'Syn', 'Red', 'Unq1', 'Unq2'};
elseif ismember('PID_atoms', reqOutputs)
    reqOutputs(ismember(reqOutputs, 'PID_atoms')) = [];
    reqOutputs = [reqOutputs, {'Syn', 'Red', 'Unq1', 'Unq2'}];
elseif ismember('all', reqOutputs)
    if strcmp(opts.redundancy_measure, 'I_BROJA')
        reqOutputs = {'Syn', 'Red', 'Unq1', 'Unq2', 'Unq', 'Joint', 'Union'};
    else 
        reqOutputs = {'Syn', 'Red', 'Unq1', 'Unq2', 'Unq', 'Joint', 'Union', 'q_dist'};
    end 
end
[isMember, ~] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('PID:invalidOutput', msg);
end

nSources = length(inputs)-1;
if nSources < 2
    msg = 'At least two sources are required.';
    error('PID:NotEnoughSources', msg);
end
if ismember('q_dist', reqOutputs)
    if ~strcmp(opts.bias, 'naive') || ~strcmp(opts.redundancy_measure, 'I_BROJA') || nSources > 2
        reqOutputs(ismember(reqOutputs, 'q_dist')) = [];
        warning('q_dist has been removed from reqOutputs because opts.bias must be ''naive'', redundancy_measure must be ''I_BROJA'' and not more than 2 Sources in the input.');
    end
    if isempty(reqOutputs)
        PID_values = NaN; 
        PID_naive= NaN;  
        PID_nullDist= NaN; 
        return;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Binning, reduce dimensions if necessary                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if opts.isBinned == false
    inputs_b = binning(inputs, opts);
    opts.isBinned = true;
else
    inputs_b = inputs;
end

nTimepoints_all = 1;
Dims_tmp = size(inputs_b{1});
nTrials_comp = Dims_tmp(end);
for var = 1:nSources+1
    Dims_var = size(inputs_b{var});
    if length(Dims_var) == 2
        nTrials = Dims_var(end);
        if nTrials ~= nTrials_comp
            error('PID: Inconsistent number of trials in input');
        end
    elseif length(Dims_var) == 3
        nTrials = Dims_var(end);
        if nTrials ~= nTrials_comp
            error('PID: Inconsistent number of trials in input');
        end
        nTimepoints_all = [nTimepoints_all, Dims_var(2)];
    else
        error('PID: Invalid input size');
    end
end
nTimepoints = max(nTimepoints_all);
if nTimepoints > 1
    opts.timeseries = true;
end 
for var = 1:nSources+1
    Dims_var = size(inputs_b{var});
    if length(Dims_var) == 2 && nTimepoints > 1
        inputs_b{var} = reshape(inputs_b{var}, [Dims_var(1), 1, Dims_var(2)]);
        inputs_b{var} = repmat(inputs_b{var}, [1, nTimepoints, 1]);
    end
end

inputs_1d = inputs_b;
size_tmp = size(inputs_1d{1});
nTrials_comp = size_tmp(end);
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
nullDist_opts = opts;
nullDist_opts.computeNulldist = false;
nullDist_opts.isBinned = true;

PID_nullDist = 0;
if opts.computeNulldist
    PID_nullDist = create_nullDist(inputs_b, reqOutputs, @PID, nullDist_opts);
end

corr = opts.bias;
corefunc = @PID;
if ~strcmp(corr, 'naive')
    [PID_values, PID_naive, PID_shuff_all] = correction(inputs_b, reqOutputs, corr, corefunc, opts);
    if ~opts.computeNulldist
        PID_nullDist = PID_shuff_all;
    end
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

[p_distr] = prob_estimator({inputs_1d{end}, inputs_1d{1:end-1}}, {'P(A,B,C)'}, opts);
for t = 1:nTimepoints
    if strcmp(opts.redundancy_measure,'I_BROJA')
        [PID_terms{t}, q_dist{t}] = opts.function(p_distr{t});
    else
        PID_terms{t} = opts.function(p_distr{t});
    end
end 
if opts.pid_constrained
        I1 = MI({inputs_1d{1}, inputs_1d{end}}, {'I(A;B)'}, opts);
        I2  = MI({inputs_1d{2}, inputs_1d{end}}, {'I(A;B)'}, opts);
        I12 = MI({cat(1, inputs_1d{1}, inputs_1d{2}), inputs_1d{end}}, {'I(A;B)'}, opts);
end 
PID_values = cell(1, length(reqOutputs));
for t = 1:nTimepoints
    if opts.pid_constrained
        if ~isfield(opts, 'chosen_atom')
            opts.chosen_atom = 'Syn';
        end
        if strcmp(opts.chosen_atom, 'Red')
            red = PID_terms{t}(1);
            PID_terms{t} = [red, I1{1}(t)-red, I2{1}(t)-red, I12{1}(t)-I1{1}(t)-I2{1}(t)+red];
        elseif strcmp(opts.chosen_atom, 'Unq1')
            un1 = PID_terms{t}(2);
            PID_terms{t} = [I1{1}(t)-un1, un1, I2{1}(t)-I1{1}(t)+un1, I12{1}(t)-I2{1}(t)-un1];
        elseif strcmp(opts.chosen_atom, 'Unq2')
            un2 = PID_terms{t}(3);
            PID_terms{t} = [I2{1}(t)-un2, I1{1}(t)-I2{1}(t)+un2, un2, I12{1}(t)-I1{1}(t)-un2];
        elseif strcmp(opts.chosen_atom, 'Syn')
            syn = PID_terms{t}(4);
            PID_terms{t} = [I1{1}(t)+I2{1}(t)-I12{1}(t)+syn, I12{1}(t)-I2{1}(t)-syn, I12{1}(t)-I1{1}(t)-syn , syn];
        end
    end
    for i = 1:length(reqOutputs)
        output = reqOutputs{i};
        switch output
            case 'Syn'
                PID_values{i}(1,t) = PID_terms{t}(4);
            case 'Red'
                PID_values{i}(1,t) = PID_terms{t}(1);
            case 'Unq1'
                PID_values{i}(1,t) = PID_terms{t}(2);
            case 'Unq2'
                PID_values{i}(1,t) = PID_terms{t}(3);
            case 'Unq'
                PID_values{i}(1,t) = PID_terms{t}(2)+PID_terms{t}(3);
            case 'Joint'
                PID_values{i}(1,t) = sum(PID_terms{t});
            case 'Union'
                PID_values{i}(1,t) = sum(PID_terms{t})-PID_terms{t}(4);
            case 'q_dist'
                if strcmp(opts.redundancy_measure,'I_BROJA') 
                   PID_values{i} = q_dist;                              
                else
                   PID_values{i} = NaN;
                end
        end
    end
end
PID_naive = PID_values;
if strcmp(opts.bias, 'shuffSub') && ~opts.computeNulldist
    PID_nullDist = MI_shuff_all;
else
    PID_nullDist = 0;
end
end
