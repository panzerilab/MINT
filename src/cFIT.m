function [cFIT_values, cFIT_naive, cFIT_nullDist, atom1, atom2] = cFIT(inputs, varargin)
% *function [cFIT_values, cFIT_naive, cFIT_nullDist, atom1, atom2] = cFIT(inputs, outputs, opts)*
%
% The cFIT function computes the conditioned Feature-specific Information Transfer (cFIT) values between time series data.
% Conditioned Feature-specific Information Transfer quantifies how much information about a specific feature (S)
% flows between two regions (A)(B) conditioned on region (C). 
%
% Inputs:
%   - inputs: A cell array containing the input data with N sources and one target:
%             - inputs{1}:   First data (A) with dimensions
%                            nDims X nTimepoints X nTrials
%             - inputs{2}:   Second data (B) with dimensions
%                            nDims X nTimepoints X nTrials
%             - inputs{3}:   Third data (C) with dimensions
%                            nDims X nTimepoints X nTrials
%             - inputs{4}:   Feature data (S) with dimensions
%                            nDims X (nTimepoints X) nTrials                                   
%
%   - outputs: A cell array of strings specifying which FIT values to compute:
%              - 'cFIT(A->B;S|C)' 
%              - 'cFIT(B->A;S|C)'
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
%                                  'setToZero'   : Sets NaN values to zero.
%                                  'error'       : (default) Throws an error if NaN values are detected.
%
%
% Outputs:
%   - cFIT_values: A cell array containing the computed cFIT values as specified in the outputs argument.
%   - cFIT_naive: A cell array containing the naive cFIT estimates.
%   - cFIT_nullDist: Results of the null distribution computation (0 if not performed).
%   - atom1: A cell array containing the first atom values computed during the analysis.
%   - atom2: A cell array containing the second atom values computed during the analysis.
%
% Example:
%   To compute the Feature-specific Information Transfer from two neural populations X1, X2 given X3 about S,
%   the function can be called as follows:
%   [cFIT_values, cFIT_naive, cFIT_nullDist, atom1, atom2] = cFIT({X1, X2, X3, S}, {'cFIT(A->B;S|C)', 'cFIT(B->A;S|C)'}, opts);
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
        [inputs, outputs, opts] = check_inputs('cFIT',inputs,varargin{:});
    end
else
    [inputs, outputs, opts] = check_inputs('cFIT',inputs,varargin{:});
end

nullDist_opts = opts;
nullDist_opts.computeNulldist = false;

if ~isfield(opts, 'recall')
    opts.recall = false;
end

if ~opts.recall
    nSources = length(inputs)-1;
    if nSources < 3
        msg = 'Two sources are required.';
        error('FIT:NotEnoughSources', msg);
    end
    DimsA = size(inputs{1});
    DimsB = size(inputs{2});
    DimsC = size(inputs{3});
    DimsS = size(inputs{4});
    if DimsA(end) ~= DimsB(end) || DimsA(end) ~= DimsC(end) || DimsA(end) ~= DimsS(end)
        msg = sprintf('The number of trials for A (%d), B (%d), C (%d) and S (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(3),DimsB(3),DimsC(3),DimsS(3));
        error('cFIT:InvalidInput', msg);
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

    A_delayed = zeros(DimsA(1), length(Atau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Atau) 
        index = opts.tpres{1} - Atau(tau);
        if index > 0 && index <= nTimepointsA
            A_delayed(:,tau, :) = inputs{1}(:, index, :);
        else
            msg = "The specified delay for A is not available for the timeseries.";
            error('cFIT:InvalidInput', msg);
        end
    end
    B_delayed = zeros(DimsB(1), length(Btau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Btau)
        index = opts.tpres{2} - Btau(tau);
        if index > 0 && index <= nTimepointsB
            B_delayed(:,tau, :) = inputs{2}(:, index, :);
        else
            msg = "The specified delay for B is not available for the timeseries.";
            error('cFIT:InvalidInput', msg);
        end
    end
    C_delayed = zeros(DimsC(1), length(Ctau), nTrials); %(nDims x nTaus x nTrials)
    for tau = 1:length(Ctau)
        index = opts.tpres{3} - Ctau(tau);
        if index > 0 && index <= nTimepointsC
            C_delayed(:,tau, :) = inputs{3}(:, index, :);
        else
            msg = "The specified delay for C is not available for the timeseries.";
            error('cFIT:InvalidInput', msg);
        end
    end

    A_pres = squeeze(A_delayed(:, 1, :));
    B_pres = squeeze(B_delayed(:, 1, :));
    B_past = squeeze(B_delayed(:, 2:end, :));
    A_past = squeeze(A_delayed(:, 2:end, :));
    C_past = squeeze(C_delayed);
    S = inputs{end};    
   
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
    if size(C_past,1) == nTrials
        C_past = C_past';
    end

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Step 3: Binning, reduce dimensions if necessary                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Reshape opts.bin_method
    if isscalar(opts.bin_method)
        opts.bin_method = repmat(opts.bin_method, 1, 6); 
    elseif numel(opts.bin_method) == 2
        opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 4)];
    elseif numel(opts.bin_method) == 3
        opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 2), repmat(opts.bin_method(3), 1, 2)]; 
    else
        opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 2), opts.bin_method(3), opts.bin_method(4)];  
    end

    % Reshape opts.n_bins
    if isscalar(opts.n_bins)
        opts.n_bins = repmat(opts.n_bins, 1, 6); % Repeat 5 times
    elseif numel(opts.n_bins) == 2
        opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 4)];
    elseif numel(opts.n_bins) == 3
        opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 2), repmat(opts.n_bins(3), 1, 2)];
    else
        opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 2), opts.n_bins(3), opts.n_bins(4)]; 
    end

    % Bin the data
    if ~opts.isBinned
        inputs_b = binning({A_pres,A_past,B_pres,B_past, C_past, S} ,opts);
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


possibleOutputs = {'cFIT(A->B;S|C)', 'cFIT(B->A;S|C)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('cFIT:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 4.A: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cFIT_nullDist = 0;
corr = opts.bias;
corefunc = @cFIT;
if any(opts.computeNulldist)
        nullDist_opts.recall = false;
        cFIT_nullDist = create_nullDist(inputs, outputs, @cFIT, nullDist_opts);
end
if ~strcmp(corr, 'naive')  
    opts.recall = true;
    opts.computeNulldist = false;
    [cFIT_values, cFIT_naive] = correction(inputs_1d, outputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Step 4.B: Compute Probability Distributions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
list_distributions = struct( ...
    'cFIT_AB_S_C', {{'P(B_pres,A_past,B_past,C_past,S)'}}, ...
    'cFIT_BA_S_C', {{'P(A_pres,B_past,A_past,C_past,S)'}} ...
    );

opts.multidim = true;
required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'cFIT(A->B;S|C)'
            required_distributions = [required_distributions, list_distributions.cFIT_AB_S_C{:}];
        case 'cFIT(B->A;S|C)'
            required_distributions = [required_distributions, list_distributions.cFIT_BA_S_C{:}];
    end
end

prob_dists = {};
for i = 1:length(required_distributions)
    switch required_distributions{i}
        case 'P(B_pres,A_past,B_past,C_past,S)'
            prob_dists{i} =  prob_estimator({inputs_1d{3},inputs_1d{2},inputs_1d{4},inputs_1d{5},inputs_1d{end}}, {'P(all)'}, opts);
        case 'P(A_pres,B_past,A_past,C_past,S)'
            prob_dists{i} =  prob_estimator({inputs_1d{1},inputs_1d{4},inputs_1d{2},inputs_1d{5},inputs_1d{end}}, {'P(all)'}, opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 4.C: Compute requested Output Values                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the TE values
cFIT_values = cell(1, length(outputs));
atom1 = cell(1, length(outputs));
atom2 = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'cFIT(A->B;S|C)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(B_pres,A_past,B_past,C_past,S)')};
            atoms = cfit_core(Prob_d{1}, opts);
            atoms1 = min(atoms(1,1), atoms(1,2));
            atoms2 = min(atoms(2,1), atoms(2,2));
            cFIT_values{i} = atoms1-atoms2;            
            atom1{i} = atoms(1,:);
            atom2{i} = atoms(2,:);
        case 'cFIT(B->A;S|C)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(A_pres,B_past,A_past,C_past,S)')};
            atoms = cfit_core(Prob_d{1}, opts);
            cFIT_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
    end
end
cFIT_naive = cFIT_values;
end