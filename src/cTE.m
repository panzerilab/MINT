function cTE_values = cTE(inputs, varargin)
%%% *function [TE_values, TE_naive, TE_shuff_all] = TE(inputs, outputs, varargin)*
%%%
%%% The cTE function computes conditioned Transfer Entropy (cTE) values from time series data.
%%% Conditioned Transfer entropy from a process A to another process B is the amount of uncertainty
%%% reduced in future values of B by knowing the past values of A given past values of B and C.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input time series data. Each cell represents a time series, where:
%%%             - inputs{1}: First time series (A) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%             - inputs{2}: Second time series (B) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%             - inputs{3}: Third time series (C) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%   - outputs: A cell array of strings specifying which entropies to compute.
%%%              Possible outputs include:
%%%               - 'cTE'      : Standard cond. Transfer Entropy.
%%%               - 'shcTE'    : Shuffled cond. Transfer Entropy.
%%%               - 'normcTE'  : Normalized cond. Transfer Entropy.
%%%               - 'normShcTE': Normalized Shuffled cond. Transfer Entropy.
%%%
%%%   - varargin: Optional arguments, passed as a structure. Fields may include:
%%%              - tau:            Integer that specifies the delay for the time series.
%%%                                This can be defined separately for A and B as a cell:
%%%                                opts.tau{1} = delay for A;
%%%                                opts.tau{2} = delay for B;
%%%                                opts.tau{3} = delay for C;
%%%                                If a single integer or a single cell field is provided,
%%%                                the same delay will be applied to both A and B.
%%%                                Default is 1 timepoint.
%%%
%%%              - tpres:          Specifies the present timepoint for the calculation.
%%%                                This can also be defined separately for A and B in the same manner as tau:
%%%                                opts.tpres{1} = present for A;
%%%                                opts.tpres{2} = present for B;
%%%                                opts.tpres{3} = present for C;
%%%                                If only a single value or cell is given, it applies to both.
%%%                                Default is set to the last timepoint of input A.
%%%
%%%              - singleTimepoint:Boolean (true/false) indicating whether to compute Transfer Entropy
%%%                                for a single timepoint or for the entire time series.
%%%                                If set to true, only single timepoints will be considered in the analysis.
%%%                                Default is false.
%%%
%%%              - bias:           Specifies the bias correction method to be used.
%%%                                Possible values include:
%%%                                'naive'                      :(default) - No correction applied.
%%%                                'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%%%                                'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%%%                                'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%%%                                'pt'                         :Panzeri-Treves bias correction.
%%%                                'bub'                        :
%%%                                Users can also define their own custom bias correction method
%%%                                (see help for correction.m)
%%%
%%%              - bin_method:     Cell array specifying the binning method.
%%%                                It can have one or two entries:
%%%                                If one entry is provided, it will be applied to both A, B and C.
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
%%%                                If one entry is provided, it will be used for both A, B and C.
%%%                                This integer defines how the continuous values will be
%%%                                discretized into bins for analysis.
%%%                                Default number of bins is {3}.
%%%
%%%              - suppressWarnings: Boolean (true/false) to suppress warning messages.
%%%                                  Default is false, meaning warnings will be shown.
%%% Outputs:
%%%   - cTE_values: A cell array containing the computed cTE values as specified in the outputs argument.
%%%   - cTE_naive: A cell array containing the naive cTE estimates.
%%%   - cTE_shuff_all: A value indicating the all results of the shuffling procedure (0 if not performed).
%%%
%%% Note:
%%% Input A, B and C can represent multiple neurons concatenated along the first dimension.
%%% This means that each neuron can contribute its activity data, allowing the analysis
%%% of interactions between different neurons and their influence on other time series.
%%%
%%% EXAMPLE
%%% Suppose we have two time series of groups of neurons X1 and X2, two time series of groups of neurons Y1 and Y2 and
%%% two time series of groups of neurons Z1 and Z2.
%%% (Structur data is nNeurons x nTimepoints x nTrials)
%%%
%%% We can structure our inputs as follows:
%%% Thus, the total input for A and B would be:
%%% A = cat(1, X1, X2);  % Concatenates X1 and X2 along the first dimension (neurons)
%%% B = cat(1, Y1, Y2);  % Concatenates Y1 and Y2 along the first dimension (neurons)
%%% C = cat(1, Z1, Z2);  % Concatenates Z1 and Z2 along the first dimension (neurons)
%%%
%%% To compute the Transfer Entropy from time series A (X) to B (Y) conditioned on c(Z) the function can be called as:
%%% cTE_values = cTE({A, B, C}, {'cTE'}, opts);
%%% 'opts' represents additional options you may want to include, such as
%%% specifying the delay (tau), number of bins (n_bins), and other parameters as needed.

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
        [outputs, opts] = check_inputs('cTE',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('cTE',inputs,varargin{:});
end

possibleOutputs = {'TE(A->B|C)', 'TE(B->A|C)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
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
    tPoints = min(minA_tp,minB_tp, minC_tp);

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
    C_past = C_delayed;
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
    C_past = C_delayed;
    S = inputs{end};
    A_pres = reshape(A_delayed(:, 1, :), 1, DimsA(3));  % Shape [nDims, nTrials]
    B_pres = reshape(B_delayed(:, 1, :), 1, DimsB(3));  % Shape [nDims, nTrials]
    A_past = reshape(A_delayed(:, 2:end, :), DimsA(1), size(A_delayed, 2)-1, DimsA(3));  % Shape [nDims-1, delays, nTrials]
    B_past = reshape(B_delayed(:, 2:end, :), DimsB(1), size(B_delayed, 2)-1, DimsB(3));  % Shape [nDims-1, delays, nTrials]
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
    if sizeVar(1) > 1  0 
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
    end
    if sizeVar(2) > 1
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 2);
    end
    inputs_1d{var} = squeeze(inputs_1d{var});
    inputs_1d{var} = reshape(inputs_1d{var}, 1, sizeVar(3));
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
cTE_values = cell(1, length(outputs));
cTE_naive = cell(1, length(outputs));
cTE_shuff_all = cell(1, length(outputs));
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

end



