function [cTE_values, cTE_naive, cTE_nullDist] = cTE(inputs, varargin)
% *function [cTE_values, cTE_naive, cTE_nullDist] = cTE(inputs, reqOutputs, opts)*
%
% The cTE function computes transfer entropy (cTE) between two 
% time series (A and B) conditioned on a third time series (C). The 
% function supports different configurations for the time series data and 
% allows for multiple time points, lag adjustments, and various binning 
% methods.
%
% Inputs:
%   - inputs: A cell array containing the input time series data. Each cell represents a time series, where:
%             - inputs{1}: First time series (A) with dimensions
%                          nDims X nTimepoints X nTrials
%             - inputs{2}: Second time series (B) with dimensions
%                          nDims X nTimepoints X nTrials
%             - inputs{3}: Third time series (C) with dimensions
%                          nDims X nTimepoints X nTrials
%
%   - reqOutputs: A cell array of strings specifying which cTE measures to compute.
%              Possible reqOutputs include:
%               - 'TE(A->B|C)' : Conditional Transfer Entropy from A to B given C
%               - 'TE(B->A|C)' : Conditional Transfer Entropy from B to A given C
%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - singleTimepoint:    Boolean (true/false) indicating whether to compute Transfer Entropy
%                                    for a single timepoint or the delayed timeseries
%
%              - tau:                Integer that specifies the delay for the time series.
%                                    This can be defined separately for A, B and C as a cell:
%                                    opts.tau{1} = delay for A;
%                                    opts.tau{2} = delay for B;
%                                    opts.tau{3} = delay for C;
%                                    If a single integer or a single cell field is provided,
%                                    the same delay will be applied to A, B and C (default: 1)
%
%              - tpres:              Specifies the present timepoint for the calculation.
%                                    This can also be defined separately for A, B and C in the same manner as tau:
%                                    opts.tpres{1} = present for A;
%                                    opts.tpres{2} = present for B;
%                                    opts.tpres{C} = present for C;
%                                    If only a single value or cell is given, it applies to A, B and C.
%                                    Default is length of A.
%
%              - bias:               Specifies the bias correction method to be used.
%                                    'naive'                      :(default) - No correction applied.
%                                    'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%                                    'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%                                    'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%                                    'pt'                         :Panzeri-Treves bias correction (Panzeri and Treves 1996).
%                                    'bub'                        :best upper bound(Paninsky, 2003)
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
%                                    Default number of bins is {3}
%
%              - computeNulldist:  If set to true, generates a null distribution
%                                  based on the specified inputs and core function.
%                                  When this option is enabled, the following can be specified:
%                                   - `n_samples`: The number of null samples to generate (default: 100).
%                                   - 'shuffling': Additional shuffling options to determine the variables to be 
%                                      shuffled during the computation of the null distribution (default: {'A'}).
%                                      (type 'help hShuffle' for more information).
% 
%              - suppressWarnings:  Boolean (true/false) to suppress warning messages.
%                                   Default is false, meaning warnings will be shown.
%
%              - NaN_handling:     Specifies how NaN values should be handled in the data.
%                                  Options include:
%                                  'removeTrial' : Removes trials containing NaN in any variable 
%                                                  from all input data.
%                                  'error'       : (default) Throws an error if NaN values are detected.
%
% Outputs:
%   - cTE_values: A cell array containing the computed cTE values as specified in the reqOutputs argument.
%   - cTE_naive: A cell array containing the naive cTE estimates.
%   - cTE_shuff_all: Results of the null distribution computation (0 if not performed).
%
% Note:
% Input A, B and C can represent multiple neurons concatenated along the first dimension.
% This means that each neuron can contribute its activity data, allowing the analysis
% of interactions between different neurons and their influence on other time series.
%
% EXAMPLE
% Suppose we have two time series of groups of neurons X1 and X2, two time series of groups of neurons Y1 and Y2 and
% two time series of groups of neurons Z1 and Z2.
% (Structur data is nNeurons x nTimepoints x nTrials)
%
% We can structure our inputs as follows:
% X = cat(1, X1, X2);  Concatenates X1 and X2 along the first dimension (neurons)
% Y = cat(1, Y1, Y2);  Concatenates Y1 and Y2 along the first dimension (neurons)
% Z = cat(1, Z1, Z2);  Concatenates Z1 and Z2 along the first dimension (neurons)
%
% To compute the Transfer Entropy from time series X to Y conditioned on Z the function can be called as:
% cTE_values = cTE({X, Y, Z}, {'TE(A->B|C)'}, opts);
%
% Here, 'opts' represents additional options you may want to include (see varargin options)

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
    error('cTE:notEnoughInput', msg);
end

if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            reqOutputs = varargin{1};
        end
    else
        [inputs, reqOutputs, opts] = check_inputs('cTE',inputs,varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('cTE',inputs,varargin{:});
end

possibleOutputs = {'TE(A->B|C)', 'TE(B->A|C)'};
[isMember, indices] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('cTE:invalidOutput', msg);
end

DimsA = size(inputs{1});
DimsB = size(inputs{2});
DimsC = size(inputs{3});
if DimsA(3) ~= DimsB(3) || DimsA(3) ~= DimsC(3)
    msg = sprintf('The number of trials for A (%d), B (%d) and C (%d) are not consistent. Ensure all variables have the same number of trials.', DimsA(3), DimsB(3), DimsC(3));
    error('cTE:InvalidInput', msg);
end
nTrials = DimsA(3);

nTimepointsA = DimsA(2);
nTimepointsB = DimsB(2);
nTimepointsC = DimsC(2);

if iscell(opts.tau)
    if numel(opts.tau) == 1
        opts.tau{2} = opts.tau{1};
        opts.tau{3} = opts.tau{1};
    elseif  numel(opts.tau) == 2
        opts.tau{3} = opts.tau{2};
    end
else
    opts.tau = {opts.tau};
    opts.tau{2} = opts.tau{1};
    opts.tau{3} = opts.tau{1};
end

if iscell(opts.tpres)
    if numel(opts.tpres) == 1
        opts.tpres{2} = opts.tpres{1};
        opts.tpres{3} = opts.tpres{1};
    elseif  numel(opts.tpres) == 2
        opts.tpres{3} = opts.tpres{2};
    end
else
    opts.tpres = {opts.tpres};
    opts.tpres{2} = opts.tpres{1};
    opts.tpres{3} = opts.tpres{1};
end

Atau = opts.tau{1};
Btau = opts.tau{2};
Atau = [0; Atau(:)];
Btau = [0; Btau(:)];
Ctau = opts.tau{3};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 2: Prepare Data (Shift based on tau,tpres/bin/reduce dims)        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~opts.singleTimepoint
    presA = opts.tpres{1};
    presB = opts.tpres{2};
    presC = opts.tpres{3};

    minA_tp = presA - max(Atau);
    minB_tp = presB - max(Btau);
    minC_tp = presC - max(Ctau);
    tPoints = min(minA_tp,minB_tp);

    A_delayed =  zeros(DimsA(1), length(Atau), tPoints, nTrials);
    A_delayed(:,1,:,:) = inputs{1}(:,(presA-tPoints+1):presA,:);
    for tau = 2:length(Atau)
        start_p = (presA-Atau(tau))-tPoints+1;
        end_p = (presA-Atau(tau));
        A_delayed(:,tau,:,:) = inputs{1}(:,start_p:end_p,:);
    end
    B_delayed =  zeros(DimsB(1), length(Btau), tPoints, nTrials);
    B_delayed(:,1,:,:) = inputs{1}(:,(presB-tPoints+1):presB,:);
    for tau = 2:length(Btau)
        start_p = (presB-Btau(tau))-tPoints+1;
        end_p = (presB-Btau(tau));
        B_delayed(:,tau,:,:) = inputs{1}(:,start_p:end_p,:);
    end
    C_delayed =  zeros(DimsC(1), length(Ctau), tPoints, nTrials);
    for tau = 1:length(Ctau)
        start_p = (presC-Ctau(tau))-tPoints+1;
        end_p = (presC-Ctau(tau));
        C_delayed(:,tau,:,:) = inputs{1}(:,start_p:end_p,:);
    end
    A_pres = A_delayed(:, 1, :);
    B_pres = B_delayed(:, 1, :);
    B_past = B_delayed(:, 2:end, :);
    A_past = A_delayed(:, 2:end, :);
    C_past = C_delayed(:,:,:);
    S = inputs{end};
else   
    A_delayed = zeros(DimsA(1), length(Atau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Atau)
        index = opts.tpres{1} - Atau(tau);
        if index > 0 && index <= nTimepointsA
            A_delayed(:,tau, :) = inputs{1}(:, index, :);
        else
            msg = "The specified delay for A is not available for the timeseries.";
            error('cTE:InvalidInput', msg);
        end
    end
    B_delayed = zeros(DimsB(1), length(Btau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Btau)
        index = opts.tpres{2} - Btau(tau);
        if index > 0 && index <= nTimepointsB
            B_delayed(:,tau, :) = inputs{2}(:, index, :);
        else
            msg = "The specified delay for B is not available for the timeseries.";
            error('cTE:InvalidInput', msg);
        end
    end
    C_delayed = zeros(DimsC(1), length(Ctau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Ctau)
        index = opts.tpres{3} - Ctau(tau);
        if index > 0 && index <= nTimepointsC
            C_delayed(:,tau, :) = inputs{3}(:, index, :);
        else
            msg = "The specified delay for C is not available for the timeseries.";
            error('cTE:InvalidInput', msg);
        end
    end
    A_pres = squeeze(A_delayed(:, 1, :));
    B_pres = squeeze(B_delayed(:, 1, :));
    B_past = squeeze(B_delayed(:, 2:end, :));
    A_past = squeeze(A_delayed(:, 2:end, :));
    S = inputs{end};
    A_pres = reshape(A_delayed(:, 1, :), DimsA(1), DimsA(3));  % Shape [nDims, nTrials]
    B_pres = reshape(B_delayed(:, 1, :), DimsB(1), DimsB(3));  % Shape [nDims, nTrials]

    A_past = reshape(A_delayed(:, 2:end, :), DimsA(1), size(A_delayed, 2)-1, DimsA(3));  % Shape [nDims-1, delays, nTrials]
    B_past = reshape(B_delayed(:, 2:end, :), DimsB(1), size(B_delayed, 2)-1, DimsB(3));  % Shape [nDims-1, delays, nTrials]
    C_past = C_delayed;
end  


% Reshape opts.bin_method
if numel(opts.bin_method) == 1
    opts.bin_method = repmat(opts.bin_method, 1, 5); % Repeat 5 times
elseif numel(opts.bin_method) == 2
    opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 3)]; % First 2x, second 3x
else
    opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 2), opts.bin_method(3)];  % 2x first, 2x second, 1x third
end

% Reshape opts.n_bins
if numel(opts.n_bins) == 1
    opts.n_bins = repmat(opts.n_bins, 1, 5); % Repeat 5 times
elseif numel(opts.n_bins) == 2
    opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 3)]; % First 2x, second 3x
else
    opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 2), opts.n_bins(3)];  % 2x first, 2x second, 1x third
end

% Bin the data
inputs_b = binning({A_pres,A_past,B_pres,B_past, C_past} ,opts);
opts.isBinned = true;
inputs_1d = inputs_b;
% Reduce the Dimensions if necessary
for var = 1:length(inputs_b)
    sizeVar = size(inputs_1d{var});
    if sizeVar(1) > 1
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
    end
    if length(sizeVar) > 2 && sizeVar(2) > 1
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 2);
    end
    inputs_1d{var} = squeeze(inputs_1d{var});
    inputs_1d{var} = reshape(inputs_1d{var}, 1, sizeVar(end));
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 3: Compute required Entropies for the requested Outputs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_dependencies = struct( ...
    'cTE_AB', {{'H(B_pres|B_past,C_past)', 'H(B_pres|B_past,A_past,C_past)'}}, ...
    'cTE_BA', {{'H(A_pres|A_past,C_past)', 'H(A_pres|A_past,B_past,C_past)'}} ...
    );

required_entropies = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'TE(A->B|C)'
            required_entropies = [required_entropies, entropy_dependencies.cTE_AB{:}];
        case 'TE(B->A|C)'
            required_entropies = [required_entropies, entropy_dependencies.cTE_BA{:}];
    end
end
required_entropies = unique(required_entropies);
H_values = cell(1, length(required_entropies));
H_naive = cell(1, length(required_entropies));
H_shuff_all = cell(1, length(required_entropies));

opts_entropy = opts;
opts_entropy.compute_nulldist = false;
for i = 1:length(required_entropies)
    switch required_entropies{i}
        case 'H(B_pres|B_past,C_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] =  H({inputs_1d{3}, cat(1, inputs_1d{4},inputs_1d{5})}, {'H(A|B)'}, opts_entropy);
        case 'H(B_pres|B_past,A_past,C_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{3}, cat(1, inputs_1d{4},inputs_1d{2}, inputs_1d{5})}, {'H(A|B)'}, opts_entropy);
        case 'H(A_pres|A_past,C_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] =  H({inputs_1d{1}, cat(1, inputs_1d{2},inputs_1d{5})}, {'H(A|B)'}, opts_entropy);
        case 'H(A_pres|A_past,B_past,C_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{1}, cat(1, inputs_1d{2},inputs_1d{4}, inputs_1d{5})}, {'H(A|B)'}, opts_entropy);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 4: Compute requested Output Values                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the cTE values
cTE_values = cell(1, length(reqOutputs));
cTE_naive = cell(1, length(reqOutputs));
cTE_shuff_all = cell(1, length(reqOutputs));
nOut = nargout;
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'TE(A->B|C)'
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_Bpres_BpastCpast_naive = H_values{strcmp(required_entropies, 'H(B_pres|B_past,C_past)')};
                H_Bpres_BpastApastCpast_naive = H_values{strcmp(required_entropies, 'H(B_pres|B_past,A_past,C_past)')};
                cTE_naive{i} = H_Bpres_BpastCpast_naive{1} - H_Bpres_BpastApastCpast_naive{1};
                H_Bpres_BpastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past,C_past)')});
                H_Bpres_BpastApastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past,A_past,C_past)')});
                for shuffIdx = 1:opts.shuff
                    H_Bpres_BpastCpast_shuff=  H_Bpres_BpastCpast_shuff_all(shuffIdx);
                    H_Bpres_BpastApastCpast_shuff = H_Bpres_BpastApastCpast_shuff_all(shuffIdx);
                    cTE_shuff_all{i}(shuffIdx) = H_Bpres_BpastCpast_shuff - H_Bpres_BpastApastCpast_shuff;
                end
                cTE_values{i} = cTE_naive{i} - mean(cTE_shuff_all{i});
                nOut = 1;
            else
                H_Bpres_BpastCpast = H_values{strcmp(required_entropies, 'H(B_pres|B_past,C_past)')};
                H_Bpres_BpastApastCpast = H_values{strcmp(required_entropies, 'H(B_pres|B_past,A_past,C_past)')};
                cTE_values{i} = H_Bpres_BpastCpast{1} - H_Bpres_BpastApastCpast{1};
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_Bpres_BpastCpast_naive = H_values{strcmp(required_entropies, 'H(B_pres|B_past,C_past)')};
                H_Bpres_BpastApastCpast_naive = H_values{strcmp(required_entropies, 'H(B_pres|B_past,A_past,C_past)')};
                cTE_naive{i} = H_Bpres_BpastCpast_naive{1} - H_Bpres_BpastApastCpast_naive{1};
                if nOut > 2 && strcmp(opts.bias,'shuffSub')
                    H_Bpres_BpastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past,C_past)')});
                    H_Bpres_BpastApastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past,A_past,C_past)')});
                    for shuffIdx = 1:opts.shuff
                        H_Bpres_BpastCpast_shuff=  H_Bpres_BpastCpast_shuff_all(shuffIdx);
                        H_Bpres_BpastApastCpast_shuff = H_Bpres_BpastApastCpast_shuff_all(shuffIdx);
                        cTE_shuff_all{i}(shuffIdx) = H_Bpres_BpastCpast_shuff - H_Bpres_BpastApastCpast_shuff;
                    end
                end
            end
        case 'TE(B->A|C)'
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_Apres_ApastCpast_naive = H_values{strcmp(required_entropies, 'H(A_pres|A_past,C_past)')};
                H_Apres_ApastBpastCpast_naive = H_values{strcmp(required_entropies, 'H(A_pres|A_past,B_past,C_past)')};
                cTE_naive{i} = H_Apres_ApastCpast_naive{1} - H_Apres_ApastBpastCpast_naive{1};
                H_Apres_ApastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past,C_past)')});
                H_Apres_ApastBpastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past,B_past,C_past)')});
                for shuffIdx = 1:opts.shuff
                    H_Apres_ApastCpast_shuff=  H_Apres_ApastCpast_shuff_all(shuffIdx);
                    H_Apres_ApastBpastCpast_shuff = H_Apres_ApastBpastCpast_shuff_all(shuffIdx);
                    cTE_shuff_all{i}(shuffIdx) = H_Apres_ApastCpast_shuff - H_Apres_ApastBpastCpast_shuff;
                end
                cTE_values{i} = cTE_naive{i} - mean(cTE_shuff_all{i});
                nOut = 1;
            else
                H_Apres_ApastCpast = H_values{strcmp(required_entropies, 'H(A_pres|A_past,C_past)')};
                H_Apres_ApastBpastCpast = H_values{strcmp(required_entropies, 'H(A_pres|A_past,B_past,C_past)')};
                cTE_values{i} = H_Apres_ApastCpast{1} - H_Apres_ApastBpastCpast{1};
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_Apres_ApastCpast_naive = H_values{strcmp(required_entropies, 'H(A_pres|A_past,C_past)')};
                H_Apres_ApastBpastCpast_naive = H_values{strcmp(required_entropies, 'H(A_pres|A_past,B_past,C_past)')};
                cTE_naive{i} = H_Apres_ApastCpast_naive{1} - H_Apres_ApastBpastCpast_naive{1};
                if nOut > 2 && strcmp(opts.bias,'shuffSub')
                    H_Apres_ApastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past,C_past)')});
                    H_Apres_ApastBpastCpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past,B_past,C_past)')});
                    for shuffIdx = 1:opts.shuff
                        H_Apres_ApastCpast_shuff=  H_Apres_ApastCpast_shuff_all(shuffIdx);
                        H_Apres_ApastBpastCpast_shuff = H_Apres_ApastBpastCpast_shuff_all(shuffIdx);
                        cTE_shuff_all{i}(shuffIdx) = H_Apres_ApastCpast_shuff - H_Apres_ApastBpastCpast_shuff;
                    end
                end
            end
    end
end
if opts.computeNulldist
    nullDist_opts = opts;
    nullDist_opts.computeNulldist = false;
    cTE_nullDist = create_nullDist(inputs, reqOutputs, @TE, nullDist_opts);
else
    cTE_nullDist = cTE_shuff_all;
end

end



