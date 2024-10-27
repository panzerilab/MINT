function [TE_values, TE_naive, TE_nullDist] = TE(inputs, varargin)
% *function [TE_values, TE_naive, TE_nullDist] = TE(inputs, outputs, opts)*
%
% The TE function computes Transfer Entropy (TE) values from time series data.
% Transfer entropy from a process A to another process B is the amount of uncertainty
% reduced in future values of B by knowing the past values of A given past values of B.
%
% Inputs:
%   - inputs: A cell array containing the data:
%             - inputs{1}: First input data (A) with dimensions
%                          nDims X nTimepoints X nTrials
%             - inputs{2}: Second input data (B) with dimensions
%                          nDims X nTimepoints X nTrials
%
%   - outputs: A cell array of strings specifying which transfer entropies to compute:
%              - 'TE(A->B;S)' : Transfer entropy from A to B
%              - 'TE(B->A;S)' : Transfer entropy from B to A
%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - singleTimepoint:    Boolean (true/false) indicating whether to compute Transfer Entropy
%                                    for a single timepoint or the delayed timeseries
%
%              - tau:                Integer that specifies the delay for the time series.
%                                    This can be defined separately for A and B as a cell:
%                                    opts.tau{1} = delay for A;
%                                    opts.tau{2} = delay for B;
%                                    If a single integer or a single cell field is provided,
%                                    the same delay will be applied to both A and B (default: 1)
%
%              - tpres:              Specifies the present timepoint for the calculation.
%                                    This can also be defined separately for A and B in the same manner as tau:
%                                    opts.tpres{1} = present for A;
%                                    opts.tpres{2} = present for B;
%                                    If only a single value or cell is given, it applies to both.
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
%                                      (type 'help shuffle' for more information).
% 
%              - suppressWarnings:  Boolean (true/false) to suppress warning messages.
%                                   Default is false, meaning warnings will be shown.
%
%              - NaN_handling:     Specifies how NaN values should be handled in the data.
%                                  Options include:
%                                  'removeTrial' : Removes trials containing NaN in any variable 
%                                                  from all input data.
%                                  'setToZero'   : Sets NaN values to zero.
%                                  'error'       : (default) Throws an error if NaN values are detected.
%
% Outputs:
%   - TE_values: A cell array containing the computed TE values as specified in the outputs argument.
%   - TE_naive: A cell array containing the naive TE estimates.
%   - TE_nullDist: Results of the null distribution computation (0 if not performed).
%
% Note:
% Input A and B can represent multiple neurons concatenated along the first dimension.
%
% EXAMPLE
% Suppose we have two time series of groups of neurons X1 and X2
% and two time series of groups of neurons Y1 and Y2.
% (Structur of X1, X2, Y1, Y1 is nNeurons x nTimepoints x nTrials)
%
% We can structure our inputs as follows:
% Thus, the total input for A and B would be:
% X = cat(1, X1, X2);   Concatenates X1 and X2 along the first dimension (neurons)
% Y = cat(1, Y1, Y2);   Concatenates Y1 and Y2 along the first dimension (neurons)
%
% To compute the Transfer Entropy from time series A (X) to B (Y), the function can be called as:
% TE_values = TE({X, Y}, {'TE(A->B;S)', 'TE(B->A;S)'}, opts);
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
if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            outputs = varargin{1};
        end
    else
        [inputs, outputs, opts] = check_inputs('TE',inputs,varargin{:});
    end
else
    [inputs, outputs, opts] = check_inputs('TE',inputs,varargin{:});
end


DimsA = size(inputs{1});
DimsB = size(inputs{2});
if DimsA(end) ~= DimsB(end)
    msg = sprintf('The number of trials for A (%d) and B (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(3),DimsB(3));
    error('TE:InvalidInput', msg);
end
nTrials = DimsA(3);

nTimepointsA = DimsA(2);
nTimepointsB = DimsB(2);

if iscell(opts.tau)
    if numel(opts.tau) == 1
        opts.tau{2} = opts.tau{1};
    end
else
    opts.tau = {opts.tau};
    opts.tau{2} = opts.tau{1};
end
if iscell(opts.tpres)
    if numel(opts.tpres) == 1
        opts.tpres{2} = opts.tpres{1};
    end
else
    opts.tpres = {opts.tpres};
    opts.tpres{2} = opts.tpres{1};
end
for i = 1:length(opts.tpres)
    if (opts.tpres{i} - abs(opts.tau{i})) <= 0
        error('tPres (%d) cannot be equal or smaller than the corresponding delay value tau (%d)', opts.tpres{i}, abs(opts.tau{i}));
    end
end

Atau = opts.tau{1};
Btau = opts.tau{2};
Atau = [0; Atau(:)];
Btau = [0; Btau(:)];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 2: Prepare Data (Shift based on tau,tpres/bin/reduce dims)        %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~opts.singleTimepoint
    presA = opts.tpres{1};
    presB = opts.tpres{2};

    minA_tp = presA - max(Atau);
    minB_tp = presB - max(Btau);
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
    A_pres = A_delayed(:, 1, :);
    B_pres = B_delayed(:, 1, :);
    B_past = B_delayed(:, 2:end, :);
    A_past = A_delayed(:, 2:end, :);
    S = inputs{end};
else   
    A_delayed = zeros(DimsA(1), length(Atau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Atau)
        index = opts.tpres{1} - Atau(tau);
        if index > 0 && index <= nTimepointsA
            A_delayed(:,tau, :) = inputs{1}(:, index, :);
        else
            msg = "The specified delay for A is not available for the timeseries.";
            error('TE:InvalidInput', msg);
        end
    end
    B_delayed = zeros(DimsB(1), length(Btau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Btau)
        index = opts.tpres{2} - Btau(tau);
        if index > 0 && index <= nTimepointsB
            B_delayed(:,tau, :) = inputs{2}(:, index, :);
        else
            msg = "The specified delay for B is not available for the timeseries.";
            error('TE:InvalidInput', msg);
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
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 3: Binning, reduce dimensions if necessary                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape opts.bin_method
if isscalar(opts.bin_method)
    opts.bin_method = repmat(opts.bin_method, 1, 4);
else
    opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 2)];
end

% Reshape opts.n_bins
if isscalar(opts.n_bins)
    opts.n_bins = repmat(opts.n_bins, 1, 4);
else
    opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 2)];
end

% Bin the data
inputs_b = binning({A_pres,A_past,B_pres,B_past} ,opts);
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

possibleOutputs = {'TE(A->B)', 'TE(B->A)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('TE:invalidOutput', msg);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 3: Compute required Entropies for the requested Outputs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% List of required entropy dependencies
entropy_dependencies = struct( ...
    'TE_AB', {{'H(B_pres|B_past)', 'H(B_pres|B_past,A_past)'}}, ...
    'TE_BA', {{'H(A_pres|A_past)', 'H(A_pres|A_past,B_past)'}} ...    
    );

required_entropies = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'TE(A->B)'
            required_entropies = [required_entropies, entropy_dependencies.TE_AB{:}];
        case 'TE(B->A)'
            required_entropies = [required_entropies, entropy_dependencies.TE_BA{:}];        
    end
end
required_entropies = unique(required_entropies);
H_values = cell(1, length(required_entropies));
H_naive = cell(1, length(required_entropies));
H_shuff_all = cell(1, length(required_entropies));

opts_entropy = opts;
opts_entropy.computeNulldist = false;
for i = 1:length(required_entropies)
    switch required_entropies{i}
        case 'H(B_pres|B_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{3}(1,:), inputs_1d{4}(1,:)}, {'H(A|B)'}, opts_entropy);
        case 'H(B_pres|B_past,A_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{3}(1,:), cat(1, inputs_1d{4}(1,:),inputs_1d{2}(1,:))}, {'H(A|B)'}, opts_entropy);                      
        case 'H(A_pres|A_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{1}(1,:), inputs_1d{2}(1,:)}, {'H(A|B)'}, opts_entropy);
        case 'H(A_pres|A_past,B_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{1}(1,:), cat(1, inputs_1d{2}(1,:),inputs_1d{4}(1,:))}, {'H(A|B)'}, opts_entropy);          
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 4: Compute requested Output Values                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the TE values
TE_values = cell(1, length(outputs));
TE_naive = cell(1, length(outputs));
TE_shuff_all = cell(1, length(outputs));
nOut = nargout;
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'TE(A->B)'
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_Bpres_Bpast_naive = H_naive{strcmp(required_entropies, 'H(B_pres|B_past)')};
                H_Bpres_Bpast_Apast_naive = H_naive{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')};
                TE_naive{i} = H_Bpres_Bpast_naive{1} - H_Bpres_Bpast_Apast_naive{1};
                H_Bpres_Bpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past)')});
                H_Bpres_Bpast_Apast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')});
                for shuffIdx = 1:opts.shuff
                    H_Bpres_Bpast_shuff = H_Bpres_Bpast_shuff_all(shuffIdx);
                    H_Bpres_Bpast_Apast_shuff = H_Bpres_Bpast_Apast_shuff_all(shuffIdx);
                    TE_shuff_all{i}(shuffIdx) = H_Bpres_Bpast_shuff - H_Bpres_Bpast_Apast_shuff;
                end
                TE_values{i} = TE_naive{i} - mean(TE_shuff_all{i});
                nOut = 1;
            else
                H_Bpres_Bpast = H_values{strcmp(required_entropies, 'H(B_pres|B_past)')};
                H_Bpres_Bpast_Apast = H_values{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')};
                TE_values{i} = H_Bpres_Bpast{1} - H_Bpres_Bpast_Apast{1};
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_Bpres_Bpast_naive = H_naive{strcmp(required_entropies, 'H(B_pres|B_past)')};
                H_Bpres_Bpast_Apast_naive = H_naive{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')};
                TE_naive{i} = H_Bpres_Bpast_naive{1} - H_Bpres_Bpast_Apast_naive{1};
                if nOut > 2 && strcmp(opts.bias,'shuffSub')
                    H_Bpres_Bpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past)')});
                    H_Bpres_Bpast_Apast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')});
                    for shuffIdx = 1:opts.shuff
                        H_Bpres_Bpast_shuff = H_Bpres_Bpast_shuff_all(shuffIdx);
                        H_Bpres_Bpast_Apast_shuff = H_Bpres_Bpast_Apast_shuff_all(shuffIdx);
                        TE_shuff_all{i}(shuffIdx) = H_Bpres_Bpast_shuff - H_Bpres_Bpast_Apast_shuff;
                    end
                end
            end

        case 'TE(B->A)'
            %____________________________________________________________________________________________%
            if strcmp(opts.bias, 'shuffSub')
                H_Apres_Apast_naive = H_naive{strcmp(required_entropies, 'H(A_pres|A_past)')};
                H_Apres_Apast_Bpast_naive = H_naive{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')};
                TE_naive{i} = H_Apres_Apast_naive{1} - H_Apres_Apast_Bpast_naive{1};
                H_Apres_Apast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past)')});
                H_Apres_Apast_Bpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')});
                for shuffIdx = 1:opts.shuff
                    H_Apres_Apast_shuff = H_Apres_Apast_shuff_all(shuffIdx);
                    H_Apres_Apast_Bpast_shuff = H_Apres_Apast_Bpast_shuff_all(shuffIdx);
                    TE_shuff_all{i}(shuffIdx) = H_Apres_Apast_shuff - H_Apres_Apast_Bpast_shuff;
                end
                TE_values{i} = TE_naive{i} - mean(TE_shuff_all{i});
                nOut = 1;
            else
                H_Apres_Apast = H_values{strcmp(required_entropies, 'H(A_pres|A_past)')};
                H_Apres_Apast_Bpast = H_values{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')};
                TE_values{i} = H_Apres_Apast{1} - H_Apres_Apast_Bpast{1};
            end
            %____________________________________________________________________________________________%
            if nOut > 1
                H_Apres_Apast_naive = H_naive{strcmp(required_entropies, 'H(A_pres|A_past)')};
                H_Apres_Apast_Bpast_naive = H_naive{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')};
                TE_naive{i} = H_Apres_Apast_naive{1} - H_Apres_Apast_Bpast_naive{1};
                if nOut > 2 && strcmp(opts.bias,'shuffSub')
                    H_Apres_Apast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past)')});
                    H_Apres_Apast_Bpast_shuff_all = cell2mat(H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')});
                    for shuffIdx = 1:opts.shuff
                        H_Apres_Apast_shuff = H_Apres_Apast_shuff_all(shuffIdx);
                        H_Apres_Apast_Bpast_shuff = H_Apres_Apast_Bpast_shuff_all(shuffIdx);
                        TE_shuff_all{i}(shuffIdx) = H_Apres_Apast_shuff - H_Apres_Apast_Bpast_shuff;
                    end
                end
            end        
    end
end

if opts.computeNulldist
    nullDist_opts = opts;
    nullDist_opts.computeNulldist = false;
    TE_nullDist = create_nullDist(inputs, outputs, @TE, nullDist_opts);
else
    TE_nullDist = TE_shuff_all;
end
end

