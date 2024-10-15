function [TE_values, TE_naive, TE_nullDist] = TE(inputs, varargin)
%%% *function [TE_values, TE_naive, TE_nullDist] = TE(inputs, outputs, varargin)*
%%%
%%% The TE function computes Transfer Entropy (TE) values from time series data.
%%% Transfer entropy from a process A to another process B is the amount of uncertainty
%%% reduced in future values of B by knowing the past values of A given past values of B.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input time series data. Each cell represents a time series, where:
%%%             - inputs{1}: First time series (A) with dimensions
%%%                          nDims X nTimepoints X nTrials
%%%             - inputs{2}: Second time series (B) with dimensions
%%%                          nDims X nTimepoints X nTrials)
%%%   - outputs: A cell array of strings specifying which entropies to compute.
%%%              Possible outputs include:
%%%               - 'TE'      : Standard Transfer Entropy.
%%%               - 'shTE'    : Shuffled Transfer Entropy.
%%%               - 'normTE'  : Normalized Transfer Entropy.
%%%               - 'normShTE': Normalized Shuffled Transfer Entropy.
%%%
%%%   - varargin: Optional arguments, passed as a structure. Fields may include:
%%%              - tau:            Integer that specifies the delay for the time series.
%%%                                This can be defined separately for A and B as a cell:
%%%                                opts.tau{1} = delay for A;
%%%                                opts.tau{2} = delay for B;
%%%                                If a single integer or a single cell field is provided,
%%%                                the same delay will be applied to both A and B.
%%%
%%%                                Default is 1 timepoint.
%%%              - tpres:          Specifies the present timepoint for the calculation.
%%%                                This can also be defined separately for A and B in the same manner as tau:
%%%                                opts.tpres{1} = present for A;
%%%                                opts.tpres{2} = present for B;
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
%%%   - TE_values: A cell array containing the computed TE values as specified in the outputs argument.
%%%   - TE_naive: A cell array containing the naive TE estimates.
%%%   - TE_nullDist: A value indicating the all results of the shuffling procedure (0 if not performed).
%%%
%%% Note:
%%% Input A and B can represent multiple neurons concatenated along the first dimension.
%%% This means that each neuron can contribute its activity data, allowing the analysis
%%% of interactions between different neurons and their influence on other time series.
%%%
%%% EXAMPLE
%%% Suppose we have two time series of groups of neurons X1 and X2
%%% and two time series of groups of neurons Y1 and Y2.
%%% (Structur of X1, X2, Y1, Y1 is nNeurons x nTimepoints x nTrials)
%%%
%%% We can structure our inputs as follows:
%%% Thus, the total input for A and B would be:
%%% A = cat(1, X1, X2);  % Concatenates X1 and X2 along the first dimension (neurons)
%%% B = cat(1, Y1, Y2);  % Concatenates Y1 and Y2 along the first dimension (neurons)
%%%
%%% To compute the Transfer Entropy from time series A (X) to B (Y), the function can be called as:
%%% TE_values = TE({A, B}, {'TE'}, opts);
%%%
%%% Here, 'opts' represents additional options you may want to include, such as
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
        [outputs, opts] = check_inputs('TE',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('TE',inputs,varargin{:});
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
    A_pres = reshape(A_delayed(:, 1, :), 1, DimsA(3));  % Shape [nDims, nTrials]
    B_pres = reshape(B_delayed(:, 1, :), 1, DimsB(3));  % Shape [nDims, nTrials]

    A_past = reshape(A_delayed(:, 2:end, :), DimsA(1), size(A_delayed, 2)-1, DimsA(3));  % Shape [nDims-1, delays, nTrials]
    B_past = reshape(B_delayed(:, 2:end, :), DimsB(1), size(B_delayed, 2)-1, DimsB(3));  % Shape [nDims-1, delays, nTrials]
end  

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 3: Binning, reduce dimensions if necessary                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Reshape opts.bin_method
if numel(opts.bin_method) == 1
    opts.bin_method = repmat(opts.bin_method, 1, 4);
else
    opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 2)];
end

% Reshape opts.n_bins
if numel(opts.n_bins) == 1
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
    if sizeVar(2) > 1
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 2);
    end
    inputs_1d{var} = squeeze(inputs_1d{var});
    inputs_1d{var} = reshape(inputs_1d{var}, 1, sizeVar(3));
end

possibleOutputs = {'TE(A->B)', 'TE(B->A)','shTE(A->B)','shTE(B->A)','normTE(A->B)','normTE(B->A)','normShTE(A->B)','normShTE(B->A)'};
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
    'TE_BA', {{'H(A_pres|A_past)', 'H(A_pres|A_past,B_past)'}}, ...
    'shTE_AB', {{'H(A_past|B_past)','H(A_past|B_past,B_pres)', 'Hsh(A_past|B_past,B_pres)'}}, ...
    'shTE_BA', {{'H(B_past|A_past)','H(B_past|A_past,A_pres)', 'Hsh(B_past|A_past,A_pres)'}}, ...
    'normTE_AB', {{'H(B_pres|B_past)', 'H(B_pres|B_past,A_past)', 'H(B_pres,B_past)', 'H(B_pres)'}}, ...
    'normTE_BA', {{'H(A_pres|A_past)', 'H(A_pres|A_past,B_past)', 'H(A_pres,A_past)', 'H(A_pres)'}}, ...
    'normShTE_AB', {{'H(A_past|B_past)','H(A_past|B_past,B_pres)', 'Hsh(A_past|B_past,B_pres)','H(B_pres,B_past)', 'H(B_pres)'}}, ...
    'normShTE_BA', {{'H(B_past|A_past)','H(B_past|A_past,A_pres)', 'Hsh(B_past|A_past,A_pres)','H(A_pres,A_past)', 'H(A_pres)'}} ...
    );

required_entropies = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'TE(A->B)'
            required_entropies = [required_entropies, entropy_dependencies.TE_AB{:}];
        case 'TE(B->A)'
            required_entropies = [required_entropies, entropy_dependencies.TE_BA{:}];
        case 'shTE(A->B)'
            required_entropies = [required_entropies, entropy_dependencies.shTE_AB{:}];
        case 'shTE(B->A)'
            required_entropies = [required_entropies, entropy_dependencies.shTE_BA{:}];
        case 'normTE(A->B)'
            required_entropies = [required_entropies, entropy_dependencies.normTE_AB{:}];
        case 'normTE(B->A)'
            required_entropies = [required_entropies, entropy_dependencies.normTE_BA{:}];
        case 'normShTE(A->B)'
            required_entropies = [required_entropies, entropy_dependencies.normShTE_AB{:}];
        case 'normShTE(B->A)'
            required_entropies = [required_entropies, entropy_dependencies.normShTE_BA{:}];
    end
end
required_entropies = unique(required_entropies);
H_values = cell(1, length(required_entropies));
H_naive = cell(1, length(required_entropies));
H_shuff_all = cell(1, length(required_entropies));

opts_entropy = opts;
opts_entropy.compute_nulldist = false;
%     inputs_b = binning({A_pres,A_past,B_pres,B_past} ,opts);
for i = 1:length(required_entropies)
    switch required_entropies{i}
        case 'H(B_pres|B_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{3}, inputs_1d{4}}, {'H(A|B)'}, opts_entropy);
        case 'H(B_pres|B_past,A_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{3}, cat(1, inputs_1d{4},inputs_1d{2})}, {'H(A|B)'}, opts_entropy);
        case 'H(A_past|B_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{2},inputs_1d{4}}, {'H(A|B)'}, opts_entropy);
        case 'H(A_past|B_past, B_pres)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{2}, cat(1, inputs_1d{4}, inputs_1d{3})}, {'H(A|B)'}, opts_entropy);
        case 'Hsh(A_past|B_past,B_pres)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{2}, cat(1,inputs_1d{4},inputs_1d{3})}, {'Hsh(A|B)'}, opts_entropy);
        case 'H(B_pres,B_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({cat(1,inputs_1d{3},inputs_1d{4})}, {'H(A)'}, opts_entropy);
        case 'H(B_pres)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{3}}, {'H(A)'}, opts_entropy);
        case 'H(A_pres|A_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{1}, inputs_1d{2}}, {'H(A|B)'}, opts_entropy);
        case 'H(A_pres|A_past,B_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{1}, cat(1, inputs_1d{2},inputs_1d{4})}, {'H(A|B)'}, opts_entropy);
        case 'H(B_past|A_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{4},inputs_1d{2}}, {'H(A|B)'}, opts_entropy);
        case 'H(B_past|A_past, A_pres)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{4}, cat(1, inputs_1d{2}, inputs_1d{1})}, {'H(A|B)'}, opts_entropy);
        case 'Hsh(B_past|A_past,A_pres)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{4}, cat(1,inputs_1d{2},inputs_1d{1})}, {'Hsh(A|B)'}, opts_entropy);
        case 'H(A_pres,A_past)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({cat(1,inputs_1d{1},inputs_1d{2})}, {'H(A)'}, opts_entropy);
        case 'H(A_pres)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs_1d{1}}, {'H(A)'}, opts_entropy);
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
        % case 'shTE(A->B)'
        %     H_Apast_Bpast = H_values{strcmp(required_entropies, 'H(A_past|B_past)')};
        %     H_Apast_Bpast_Bpres = H_values{strcmp(required_entropies, 'H(A_past|B_past, B_pres)')};
        %     Hsh_Apast_Bpast_Bpres = H_values{strcmp(required_entropies, 'Hsh(A_past|B_past,B_pres)')};
        %     TE_values{i} = H_Apast_Bpast{1} - H_Apast_Bpast_Bpres{1} + Hsh_Apast_Bpast_Bpres{1};
        %     if nOut > 1
        %         H_Apast_Bpast = H_naive{strcmp(required_entropies, 'H(A_past|B_past)')};
        %         H_Apast_Bpast_Bpres = H_naive{strcmp(required_entropies, 'H(A_past|B_past, B_pres)')};
        %         Hsh_Apast_Bpast_Bpres = H_naive{strcmp(required_entropies, 'Hsh(A_past|B_past,B_pres)')};
        %         TE_naive{i} = H_Apast_Bpast{1} -  H_Apast_Bpast_Bpres{1} + Hsh_Apast_Bpast_Bpres{1};
        %         if nOut > 2 && strcmp(opts.bias,'shuffSub')
        %             for shuffIdx = 1:opts.shuff
        %                 H_Apast_Bpast = H_shuff_all{strcmp(required_entropies, 'H(A_past|B_past)')}(shuffIdx);
        %                 H_Apast_Bpast_Bpres = H_shuff_all{strcmp(required_entropies, 'H(A_past|B_past, B_pres)')}(shuffIdx);
        %                 Hsh_Apast_Bpast_Bpres = H_shuff_all{strcmp(required_entropies, 'Hsh(A_past|B_past,B_pres)')}(shuffIdx);
        %                 TE_shuff_all{i}(shuffIdx) = H_Apast_Bpast{1} -  H_Apast_Bpast_Bpres{1} + Hsh_Apast_Bpast_Bpres{1};
        %             end
        %         end
        %     end
        % case 'shTE(B->A)'
        %     H_Bpast_Apast = H_values{strcmp(required_entropies, 'H(B_past|A_past)')};
        %     H_Bpast_Apast_Apres = H_values{strcmp(required_entropies, 'H(B_past|A_past, A_pres)')};
        %     Hsh_Bpast_Apast_Apres = H_values{strcmp(required_entropies, 'Hsh(B_past|A_past,A_pres)')};
        %     TE_values{i} = H_Bpast_Apast{1} - H_Bpast_Apast_Apres{1} + Hsh_Bpast_Apast_Apres{1};
        %     if nOut > 1
        %         H_Bpast_Apast = H_naive{strcmp(required_entropies, 'H(B_past|A_past)')};
        %         H_Bpast_Apast_Apres = H_naive{strcmp(required_entropies, 'H(B_past|A_past, A_pres)')};
        %         Hsh_Bpast_Apast_Apres = H_naive{strcmp(required_entropies, 'Hsh(B_past|A_past,A_pres)')};
        %         TE_naive{i} = H_Bpast_Apast{1} -  H_Bpast_Apast_Apres{1} + Hsh_Bpast_Apast_Apres{1};
        %         if nOut > 2 && strcmp(opts.bias,'shuffSub')
        %             for shuffIdx = 1:opts.shuff
        %                 H_Bpast_Apast = H_shuff_all{strcmp(required_entropies, 'H(B_past|A_past)')}(shuffIdx);
        %                 H_Bpast_Apast_Apres = H_shuff_all{strcmp(required_entropies, 'H(B_past|A_past, A_pres)')}(shuffIdx);
        %                 Hsh_Bpast_Apast_Apres = H_shuff_all{strcmp(required_entropies, 'Hsh(B_past|A_past,A_pres)')}(shuffIdx);
        %                 TE_shuff_all{i}(shuffIdx) = H_Bpast_Apast{1} -  H_Bpast_Apast_Apres{1} + Hsh_Bpast_Apast_Apres{1};
        %             end
        %         end
        %     end
        % case 'normTE(A->B)'
        %     H_Bpres_Bpast = H_values{strcmp(required_entropies, 'H(B_pres|B_past)')};
        %     H_Bpres_Bpast_Apast = H_values{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')};
        %     H_BpresBpast = H_values{strcmp(required_entropies, 'H(B_pres,B_past)')};
        %     H_Bpres = H_values{strcmp(required_entropies, 'H(B_pres)')};
        %     TE_values{i} = (H_Bpres_Bpast{1} - H_Bpres_Bpast_Apast{1})/(H_BpresBpast{1}-H_Bpres{1});
        %     if nOut > 1
        %         H_Bpres_Bpast = H_naive{strcmp(required_entropies, 'H(B_pres|B_past)')};
        %         H_Bpres_Bpast_Apast = H_naive{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')};
        %         H_BpresBpast = H_naive{strcmp(required_entropies, 'H(B_pres,B_past)')};
        %         H_Bpres = H_naives{strcmp(required_entropies, 'H(B_pres)')};
        %         TE_naive{i} = (H_Bpres_Bpast{1} - H_Bpres_Bpast_Apast{1})/(H_BpresBpast{1}-H_Bpres{1});
        %         if nOut > 2 && strcmp(opts.bias,'shuffSub')
        %             for shuffIdx = 1:opts.shuff
        %                 H_Bpres_Bpast = H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past)')}(shuffIdx);
        %                 H_Bpres_Bpast_Apast = H_shuff_all{strcmp(required_entropies, 'H(B_pres|B_past,A_past)')}(shuffIdx);
        %                 H_BpresBpast = H_shuff_all{strcmp(required_entropies, 'H(B_pres,B_past)')}(shuffIdx);
        %                 H_Bpres = H_shuff_all{strcmp(required_entropies, 'H(B_pres)')}(shuffIdx);
        %                 TE_shuff_all{i}(shuffIdx) = (H_Bpres_Bpast{1} - H_Bpres_Bpast_Apast{1})/(H_BpresBpast{1}-H_Bpres{1});
        %             end
        %         end
        %     end
        % case 'normTE(B->A)'
        %     H_Apres_Apast = H_values{strcmp(required_entropies, 'H(A_pres|A_past)')};
        %     H_Apres_Apast_Bpast = H_values{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')};
        %     H_ApresApast = H_values{strcmp(required_entropies, 'H(A_pres,A_past)')};
        %     H_Apres = H_values{strcmp(required_entropies, 'H(A_pres)')};
        %     TE_values{i} = (H_Apres_Apast{1} - H_Apres_Apast_Bpast{1})/(H_ApresApast{1}-H_Apres{1});
        %     if nOut > 1
        %         H_Apres_Apast = H_naive{strcmp(required_entropies, 'H(A_pres|A_past)')};
        %         H_Apres_Apast_Bpast = H_naive{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')};
        %         H_ApresApast = H_naive{strcmp(required_entropies, 'H(A_pres,A_past)')};
        %         H_Apres = H_naives{strcmp(required_entropies, 'H(A_pres)')};
        %         TE_naive{i} = (H_Apres_Apast{1} - H_Apres_Apast_Bpast{1})/(H_ApresApast{1}-H_Apres{1});
        %         if nOut > 2 && opts.shuff > 0
        %             for shuffIdx = 1:opts.shuff
        %                 H_Apres_Apast = H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past)')}(shuffIdx);
        %                 H_Apres_Apast_Bpast = H_shuff_all{strcmp(required_entropies, 'H(A_pres|A_past,B_past)')}(shuffIdx);
        %                 H_ApresApast = H_shuff_all{strcmp(required_entropies, 'H(A_pres,A_past)')}(shuffIdx);
        %                 H_Apres = H_shuff_all{strcmp(required_entropies, 'H(A_pres)')}(shuffIdx);
        %                 TE_shuff_all{i}(shuffIdx) = (H_Apres_Apast{1} - H_Apres_Apast_Bpast{1})/(H_ApresApast{1}-H_Apres{1});
        %             end
        %         end
        %     end
        % case 'normShTE(A->B)'
        %     H_Apast_Bpast = H_values{strcmp(required_entropies, 'H(A_past|B_past)')};
        %     H_Apast_Bpast_Bpres = H_values{strcmp(required_entropies, 'H(A_past|B_past, B_pres)')};
        %     Hsh_Apast_Bpast_Bpres = H_values{strcmp(required_entropies, 'Hsh(A_past|B_past,B_pres)')};
        %     H_BpresBpast = H_values{strcmp(required_entropies, 'H(B_pres,B_past)')};
        %     H_Bpres = H_values{strcmp(required_entropies, 'H(B_pres)')};
        %     TE_values{i} = (H_Apast_Bpast{1} -  H_Apast_Bpast_Bpres{1} + Hsh_Apast_Bpast_Bpres{1})/(H_BpresBpast-H_Bpres{1});
        %     if nOut > 1
        %         H_Apast_Bpast = H_naive{strcmp(required_entropies, 'H(A_past|B_past)')};
        %         H_Apast_Bpast_Bpres = H_naive{strcmp(required_entropies, 'H(A_past|B_past, B_pres)')};
        %         Hsh_Apast_Bpast_Bpres = H_naive{strcmp(required_entropies, 'Hsh(A_past|B_past,B_pres)')};
        %         H_BpresBpast = H_naive{strcmp(required_entropies, 'H(B_pres,B_past)')};
        %         H_Bpres = H_naive{strcmp(required_entropies, 'H(B_pres)')};
        %         TE_naive{i} = (H_Apast_Bpast{1} -  H_Apast_Bpast_Bpres{1} + Hsh_Apast_Bpast_Bpres{1})/(H_BpresBpast-H_Bpres{1});
        %         if nOut > 2 && strcmp(opts.bias,'shuffSub')
        %             for shuffIdx = 1:opts.shuff
        %                 H_Apast_Bpast = H_shuff_all{strcmp(required_entropies, 'H(A_past|B_past)')}(shuffIdx);
        %                 H_Apast_Bpast_Bpres = H_shuff_all{strcmp(required_entropies, 'H(A_past|B_past, B_pres)')}(shuffIdx);
        %                 Hsh_Apast_Bpast_Bpres = H_shuff_all{strcmp(required_entropies, 'Hsh(A_past|B_past,B_pres)')}(shuffIdx);
        %                 H_BpresBpast = H_shuff_all{strcmp(required_entropies, 'H(B_pres,B_past)')}(shuffIdx);
        %                 H_Bpres = H_shuff_all{strcmp(required_entropies, 'H(B_pres)')}(shuffIdx);
        %                 TE_shuff_all{i} = (H_Apast_Bpast{1} -  H_Apast_Bpast_Bpres{1} + Hsh_Apast_Bpast_Bpres{1})/(H_BpresBpast-H_Bpres{1});
        %             end
        %         end
        %     end
        % case 'normShTE(B->A)'
        %     H_Bpast_Apast = H_values{strcmp(required_entropies, 'H(B_past|A_past)')};
        %     H_Bpast_Apast_Apres = H_values{strcmp(required_entropies, 'H(B_past|A_past, A_pres)')};
        %     Hsh_Bpast_Apast_Apres = H_values{strcmp(required_entropies, 'Hsh(B_past|A_past,A_pres)')};
        %     H_ApresApast = H_values{strcmp(required_entropies, 'H(A_pres,A_past)')};
        %     H_Apres = H_values{strcmp(required_entropies, 'H(A_pres)')};
        %     TE_values{i} = (H_Bpast_Apast{1} -  H_Bpast_Apast_Apres{1} + Hsh_Bpast_Apast_Apres{1})/(H_ApresApast-H_Apres{1});
        %     if nOut > 1
        %         H_Bpast_Apast = H_naive{strcmp(required_entropies, 'H(B_past|A_past)')};
        %         H_Bpast_Apast_Apres = H_naive{strcmp(required_entropies, 'H(B_past|A_past, A_pres)')};
        %         Hsh_Bpast_Apast_Apres = H_naive{strcmp(required_entropies, 'Hsh(B_past|A_past,A_pres)')};
        %         H_ApresApast = H_naive{strcmp(required_entropies, 'H(A_pres,A_past)')};
        %         H_Apres = H_naive{strcmp(required_entropies, 'H(A_pres)')};
        %         TE_naive{i} = (H_Bpast_Apast{1} -  H_Bpast_Apast_Apres{1} + Hsh_Bpast_Apast_Apres{1})/(H_ApresApast-H_Apres{1});
        %         if nOut > 2 && strcmp(opts.bias,'shuffSub')
        %             for shuffIdx = 1:opts.shuff
        %                 H_Bpast_Apast = H_shuff_all{strcmp(required_entropies, 'H(B_past|A_past)')}(shuffIdx);
        %                 H_Bpast_Apast_Apres = H_shuff_all{strcmp(required_entropies, 'H(B_past|A_past, A_pres)')}(shuffIdx);
        %                 Hsh_Bpast_Apast_Apres = H_shuff_all{strcmp(required_entropies, 'Hsh(B_past|A_past,A_pres)')}(shuffIdx);
        %                 H_ApresApast = H_shuff_all{strcmp(required_entropies, 'H(A_pres,A_past)')}(shuffIdx);
        %                 H_Apres = H_shuff_all{strcmp(required_entropies, 'H(A_pres)')}(shuffIdx);
        %                 TE_shuff_all{i}(shuffIdx) = (H_Bpast_Apast{1} -  H_Bpast_Apast_Apres{1} + Hsh_Bpast_Apast_Apres{1})/(H_ApresApast-H_Apres{1});
        %             end
        %         end
            % end
    end
end

if opts.computeNulldist
    nullDist_opts = opts;
    nullDist_opts.computeNulldist = false;
    TE_nullDist = create_NullDistribution(inputs, outputs, @TE, nullDist_opts);
else
    TE_nullDist = TE_shuff_all;
end
end

