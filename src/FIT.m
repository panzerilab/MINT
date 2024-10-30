function [FIT_values, FIT_naive, FIT_nullDist, atom1, atom2] = FIT(inputs, varargin)
% *function [FIT_values, FIT_naive, FIT_nullDist, atom1, atom2] = FIT(inputs, reqOutputs, opts)*
%
% The FIT function computes the Feature-specific Information Transfer (FIT) values between time series data.
% Feature-specific Information Transfer quantifies how much information about a specific feature (S)
% flows between two regions (A)(B). 
%
% Inputs:
%   - inputs: A cell array containing the input data with N sources and one target:
%             - inputs{1}:   First data (A) with dimensions
%                            nDims X nTimepoints X nTrials
%             - inputs{2}:   Second data (B) with dimensions
%                            nDims X nTimepoints X nTrials
%             - inputs{3}:   Feature data (S) with dimensions
%                            nDims X (nTimepoints X) nTrials  
%
%   - reqOutputs: A cell array of strings specifying which FIT values to compute.
%              - 'FIT(A->B;S)' 
%              - 'FIT(B->A;S)'
%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - bias:               Specifies the bias correction method to be used.
%                                    'naive'                      :(default) - No correction applied.
%                                    'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%                                    'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%                                    'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%                                    Users can also define their own custom bias correction method
%                                    (type 'help correction' for more information)
%
%              - tau:                Cell array specifying the delays for the analysis.
%                                    Possible values include scalar or vector values indicating
%                                    the time lags to consider for A and B.
%
%              - tpres:              Cell array specifying the time points to be used as 'present'
%                                    It can have one or two entries to apply to A and B, respectively.
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
%   - FIT_values: A cell array containing the computed FIT values as specified in the reqOutputs argument.
%   - FIT_naive: A cell array containing the naive FIT estimates.
%   - FIT_nullDist: Results of the null distribution computation (0 if not performed).
%   - atom1: A cell array containing the first atom values computed during the analysis.
%   - atom2: A cell array containing the second atom values computed during the analysis.
%
% Example:
%   To compute the Feature-specific Information Transfer from two neural populations X1 and X2 about S,
%   the function can be called as follows:
%   inputs = {X1, X2, S};
%   FIT_values = FIT({X1, X2, S}s, {'FIT(A->B;S)', 'FIT(B->A;S)'}, opts);
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
            reqOutputs = varargin{1};
        end
    else
        [inputs, reqOutputs, opts] = check_inputs('FIT',inputs,varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('FIT',inputs,varargin{:});
end

nullDist_opts = opts;
nullDist_opts.computeNulldist = false;

if ~isfield(opts, 'recall')
    opts.recall = false;
end

if ~opts.recall
    nSources = length(inputs)-1;
    if nSources < 2
        msg = 'Two sources are required.';
        error('FIT:NotEnoughSources', msg);
    end
    DimsA = size(inputs{1});
    DimsB = size(inputs{2});
    DimsS = size(inputs{3});
    if DimsA(end) ~= DimsB(end) || DimsA(end) ~= DimsS(end)
        msg = sprintf('The number of trials for A (%d), B (%d) and S (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(3),DimsB(3),DimsS(3));
        error('FIT:InvalidInput', msg);
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

    Atau = opts.tau{1};
    Btau = opts.tau{2};
    Atau = [0; Atau(:)];
    Btau = [0; Btau(:)];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        Step 2: Prepare Data (Shift based on tau,tpres/bin/reduce dims)        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A_delayed = zeros(DimsA(1), length(Atau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Atau) 
        index = opts.tpres{1} - Atau(tau);
        if index > 0 && index <= nTimepointsA
            A_delayed(:,tau, :) = inputs{1}(:, index, :);
        else
            msg = "The specified delay for A is not available for the timeseries.";
            error('FIT:InvalidInput', msg);
        end
    end
    B_delayed = zeros(DimsB(1), length(Btau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Btau)
        index = opts.tpres{2} - Btau(tau);
        if index > 0 && index <= nTimepointsB
            B_delayed(:,tau, :) = inputs{2}(:, index, :);
        else
            msg = "The specified delay for B is not available for the timeseries.";
            error('FIT:InvalidInput', msg);
        end
    end

    A_pres = squeeze(A_delayed(:, 1, :));
    B_pres = squeeze(B_delayed(:, 1, :));
    B_past = squeeze(B_delayed(:, 2:end, :));
    A_past = squeeze(A_delayed(:, 2:end, :));
    S = inputs{end};

    % A_pres = reshape(A_delayed(:, 1, :), DimsA(1), DimsA(3));  % Shape [nDims, nTrials]
    % B_pres = reshape(B_delayed(:, 1, :), DimsB(1), DimsB(3));  % Shape [nDims, nTrials]
    % A_past = squeeze(reshape(A_delayed(:, 2:end, :), DimsA(1), size(A_delayed, 2)-1, DimsA(3)));  % Shape [nDims-1, delays, nTrials]
    % B_past = squeeze(reshape(B_delayed(:, 2:end, :), DimsB(1), size(B_delayed, 2)-1, DimsB(3)));  % Shape [nDims-1, delays, nTrials]
   
    if size(A_past,1) == nTrials
        A_past = A_past';
    end
    if size(A_pres,1) == nTrials
        A_pres = A_pres';
    end
    if size(B_past,1) == nTrials
        B_past = B_past';
    end
    if size(B_pres,1) == nTrials
        B_pres = B_pres';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Step 3: Binning, reduce dimensions if necessary                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

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
    if ~opts.isBinned
        inputs_b = binning({A_pres,A_past,B_pres,B_past, S} ,opts);
        opts.isBinned = true;
    end 
    inputs_1d = inputs_b;
    % Reduce the Dimensions if necessary
    for var = 1:length(inputs_b)
        sizeVar = size(inputs_1d{var});
        if sizeVar(1) > 1
            inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
        end
    end
     opts.recall = true;
else
    inputs_1d = inputs;
end


possibleOutputs = {'FIT(A->B;S)', 'FIT(B->A;S)'};
[isMember, indices] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('FIT:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 4.A: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIT_nullDist = 0;
corr = opts.bias;
corefunc = @FIT;
if any(opts.computeNulldist)
        nullDist_opts.recall = false;
        FIT_nullDist = create_nullDist(inputs, reqOutputs, @FIT, nullDist_opts);
end
if ~strcmp(corr, 'naive')  
    opts.recall = true;
    opts.computeNulldist = false;
    [FIT_values, FIT_naive] = correction(inputs_1d, reqOutputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Step 4.B: Compute Probability Distributions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_distributions = struct( ...
    'FIT_AB_S', {{'P(B_pres,A_past,B_past,S)'}}, ...
    'FIT_BA_S', {{'P(A_pres,B_Past,A_past,S)'}} ...
    );

opts.multidim = true;
required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'FIT(A->B;S)'
            required_distributions = [required_distributions, entropy_distributions.FIT_AB_S{:}];
        case 'FIT(B->A;S)'
            required_distributions = [required_distributions, entropy_distributions.FIT_BA_S{:}];
    end
end

prob_dists = {};
for i = 1:length(required_distributions)
    switch required_distributions{i}
        case 'P(B_pres,A_past,B_past,S)'
            prob_dists{i} =  prob_estimator({inputs_1d{3},inputs_1d{2},inputs_1d{4},inputs_1d{end}}, {'P(all)'}, opts);
        case 'P(A_pres,B_Past,A_past,S)'
            prob_dists{i} =  prob_estimator({inputs_1d{1},inputs_1d{4},inputs_1d{2},inputs_1d{end}}, {'P(all)'}, opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 4.C: Compute requested Output Values                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the TE values
FIT_values = cell(1, length(reqOutputs));
atom1 = cell(1, length(reqOutputs));
atom2 = cell(1, length(reqOutputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'FIT(A->B;S)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(B_pres,A_past,B_past,S)')};
            atoms = fit_core(Prob_d{1});
            FIT_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
        case 'FIT(B->A;S)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(A_pres,B_Past,A_past,S)')};
            atoms = fit_core(Prob_d{1});
            FIT_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
    end
end
FIT_naive = FIT_values;
end