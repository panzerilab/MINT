function [II_values, II_plugin, II_nullDist, atom1, atom2] = II(inputs, varargin)
% II - Intersection Information (II) and related information-theoretic quantities
%
% This function calculates the intersection information (II) based on the provided inputs, 
% reqOutputs, and optional parameters.
%
% Inputs:
%   - inputs: A cell array containing the data:
%             - inputs{1}: First input data (A) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             - inputs{2}: Second input data (B) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             - inputs{3}: Third input data (C) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             -> In cases where the input is provided as a time series, the outputs 
%                will be computed for each time point, resulting in outputs that are 
%                also represented as time series
%
%   - reqOutputs: A cell array of strings specifying which measures to compute:
%               - 'II(A,B,C)'    :information between A and B that is readout for C
%               - 'II(B,C,A)'    :information between B and C that is readout for A
%               - 'II(A,C,B)'    :information between A and C that is readout for B
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
%                                    'plugin'                      :(default) - No correction applied.
%                                    'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%                                    'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%                                    'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%                                    'ksg'                        :correction using a k-neighbors entropy estimator (Holmes and Nemenman, 2019). Only available when redundancy_measure is I_MMI
%                                    'nsb'                        :correction using the NSB algorithm (Nemenman, Bialek and van Steveninck, 2019). Only available when redundancy_measure is I_MMI
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
% Outputs:
%   - II_values: A cell array containing the computed II values as specified in the reqOutputs argument.
%   - II_plugin: A cell array containing the plugin II estimates.
%   - II_nullDist: Results of the null distribution computation (0 if not performed).
%   - atom1, atom2: The redundancy atoms between A and B about C (atom1), and between C and B about A (atom2).

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
    error('II:notEnoughInput', msg);
end

if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            reqOutputs = varargin{1};
        end
    else
        [inputs, reqOutputs, opts] = check_inputs('II',inputs,varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('II',inputs,varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Binning, reduce dimensions if necessary                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~opts.isBinned
    inputs_b = binning(inputs ,opts);
    opts.isBinned = true;
else
    inputs_b = inputs;
end

nTimepoints_all = 1;
Dims_tmp = size(inputs_b{1});
nTrials_comp = Dims_tmp(end);
for var = 1:length(inputs)
    Dims_var = size(inputs_b{var});
    if length(Dims_var) == 2
        nTrials = Dims_var(end);
        if nTrials ~= nTrials_comp
            error('II: Inconsistent number of trials in input');
        end
    elseif length(Dims_var) == 3
        nTrials = Dims_var(end);
        if nTrials ~= nTrials_comp
            error('II: Inconsistent number of trials in input');
        end
        nTimepoints_all = [nTimepoints_all, Dims_var(2)];
    else
        error('II: Invalid input size');
    end
end
nTimepoints = max(nTimepoints_all);
if nTimepoints > 1
    opts.timeseries = true;
end 
for var = 1:length(inputs)
    Dims_var = size(inputs_b{var});
    if length(Dims_var) == 2 && nTimepoints > 1
        inputs_b{var} = reshape(inputs_b{var}, [Dims_var(1), 1, Dims_var(2)]);
        inputs_b{var} = repmat(inputs_b{var}, [1, nTimepoints, 1]);
    end
end

inputs_1d = inputs_b;
size_tmp = size(inputs_1d{1});
nTrials_comp = size_tmp(end);
for var = 1:length(inputs)
    sizeVar = size(inputs_1d{var});
    nTrials = sizeVar(end);
    if nTrials ~= nTrials_comp
        msg = 'Inconsistent input size. Number of Trials must be consistent.';
        error('II:Invalid Input', msg);
    end
    if sizeVar(1) > 1
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
    end
end

possibleOutputs = {'II(A,B,C)','II(B,C,A)', 'II(A,C,B)'};
[isMember, indices] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('II:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Step 3.A: Bias correction and nullDist if requested                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II_nullDist = 0;
corr = opts.bias;
corefunc = @II;
nullDist_opts = opts;
nullDist_opts.computeNulldist = false;

if any(opts.computeNulldist)
        nullDist_opts.isBinned=true;
        II_nullDist = create_nullDist(inputs_b, reqOutputs, @II, nullDist_opts);
end
if ~strcmp(corr, 'plugin')  
    opts.computeNulldist = false;
    [II_values, II_plugin] = correction(inputs_1d, reqOutputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Step 3.B: Compute Probability Distributions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_distributions = struct( ...
    'II_ABC', {{'P(A,B,C)'}}, ...
    'II_BCA', {{'P(B,C,A)'}}, ...
    'II_ACB', {{'P(A,C,B)'}} ...
    );

switch opts.redundancy_measure
    case 'I_BROJA'
        opts.function = @pidBROJA;
    case 'I_MMI'
        opts.function = @pidimmi;
    case 'I_min'
        opts.function = @pidimin;
end

required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'II(A,B,C)'
            required_distributions = [required_distributions, entropy_distributions.II_ABC{:}];
        case 'II(B,C,A)'
            required_distributions = [required_distributions, entropy_distributions.II_BCA{:}];
        case 'II(A,C,B)'
            required_distributions = [required_distributions, entropy_distributions.II_ACB{:}];
    end
end

prob_dists = {};
for i = 1:length(required_distributions)
    switch required_distributions{i}
        case 'P(A,B,C)'
            prob_dists{i} =  prob_estimator({inputs_1d{1},inputs_1d{2},inputs_1d{3}}, {'P(A,B,C)'}, opts);
        case 'P(B,C,A)'
            prob_dists{i} =  prob_estimator({inputs_1d{2},inputs_1d{3},inputs_1d{1}}, {'P(A,B,C)'}, opts);
        case 'P(A,C,B)'
            prob_dists{i} =  prob_estimator({inputs_1d{1},inputs_1d{3},inputs_1d{2}}, {'P(A,B,C)'}, opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 3.C: Compute requested Output Values                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the TE values
II_values = cell(1, length(reqOutputs));
atom1 = cell(1, length(reqOutputs));
atom2 = cell(1, length(reqOutputs));
for t = 1:nTimepoints
    for i = 1:length(indices)
        idx = indices(i);
        switch possibleOutputs{idx}          
            case 'II(A,B,C)'
                Prob_d = prob_dists{strcmp(required_distributions, 'P(A,B,C)')};
                Prob_d = Prob_d{t, 1};
                atoms = ii_core(Prob_d, opts);
                II_values{i}(1,t) = min(atoms);
                atom1{i}(1,t) = atoms(1);
                atom2{i}(1,t) = atoms(2);
            case 'II(B,C,A)'
                Prob_d = prob_dists{strcmp(required_distributions, 'P(B,C,A)')};
                Prob_d = Prob_d{t, 1};
                atoms = ii_core(Prob_d, opts);
                II_values{i}(1,t) = min(atoms);
                atom1{i}(1,t) = atoms(1);
                atom2{i}(1,t) = atoms(2);
            case 'II(A,C,B)'
                Prob_d = prob_dists{strcmp(required_distributions, 'P(A,C,B)')};
                Prob_d = Prob_d{t, 1};
                atoms = ii_core(Prob_d, opts);
                II_values{i}(1,t) = min(atoms);
                atom1{i}(1,t) = atoms(1);
                atom2{i}(1,t) = atoms(2);
        end
    end
end 
II_plugin = II_values;
end