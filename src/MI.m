function [MI_values, MI_naive, MI_nullDist] = MI(inputs, varargin)
% *function [MI_values, MI_naive, MI_nullDist] = MI(inputs, reqOutputs, opts)*
% MI - Calculate Mutual Information (MI) and related information-theoretic quantities
%
% This function calculates mutual information (MI) and other related
% measures based on the provided inputs, requested outputs (reqOutputs), and optional parameters.
%
% Inputs:
%   - inputs: A cell array containing the data:
%             - inputs{1}: First input data (A) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             - inputs{2}: Second input data (B) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             -> In cases where the input is provided as a time series, the outputs
%                will be computed for each time point, resulting in outputs that are
%                also represented as time series
%
%   - reqOutputs: A cell array of strings specifying which entropies to compute.
%               - 'I(A;B)'    : Mutual Information I(A;B)
%               - 'Ilin(A;B)' : Linear MI Ilin(A;B)
%               - 'coI(A;B)'  : Co-information coI(A;B)
%               - 'Iss(A)'    : Sig. Sim Information Iss(A)
%               - 'Ic(A;B)'   : Sum of Ici and Icd - Ic(A;B)
%               - 'Ici(A;B)'  : Correlation Independent Information Ici(A;B)
%               - 'Icd(A;B)'  : Correlation Dependent Information Icd(A;B)
%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - bias:             Specifies the bias correction method to be used.
%                                  'naive'                      :(default) - No correction applied.
%                                  'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%                                  'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%                                  'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%                                  'pt'                         :Panzeri-Treves bias correction (Panzeri and Treves 1996).
%                                  'bub'                        :best upper bound(Paninsky, 2003)
%                                  'shuffCorr'                  :correction using stimulus-conditioned shuffling (Panzeri et al., 2007)
%                                  Users can also define their own custom bias correction method
%                                  (type 'help correction' for more information)
%
%              - bin_method:       Cell array specifying the binning method to be applied.
%                                  'none'      : (default) - No binning applied.
%                                  'eqpop'     : Equal population binning.
%                                  'eqspace'   : Equal space binning.
%                                  'userEdges' : Binning based on a specified edged.
%                                  Users can also define their own custom binning method
%                                  If one entry is provided, it will be applied to both A and B.
%                                  (type 'help binning' for more information).
%
%              - n_bins:           Specifies the number of bins to use for binning.
%                                  It can be a single integer or a cell array with one or two entries.
%                                  Default number of bins is {3}.
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
%   - MI_values: A cell array containing the computed MI values as specified in the reqOutputs argument.
%   - MI_naive: A cell array containing the naive MI estimates.
%   - MI_shuff_all: Results of the null distribution computation (0 if not performed).
%
% EXAMPLE
% Suppose we have two time series of groups of neurons X1 and X2
% and two time series of groups of neurons Y1 and Y2.
% (Structur of X1, X2, Y1, Y1 is nNeurons x nTrials)
%
% We can structure our inputs as follows:
% X = cat(1, X1, X2);  % Concatenates X1 and X2 along the first dimension (neurons)
% Y = cat(1, Y1, Y2);  % Concatenates Y1 and Y2 along the first dimension (neurons)
%
% To compute the Mutual Information and the Shuffled MI from X to Y, the function can be called as:
% MI_values = MI({X, Y}, {'I(A;B)', Ish(A;B}, opts);
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
    error('MI:notEnoughInput', msg);
end

if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            reqOutputs = varargin{1};
        end
    else
        [inputs, reqOutputs, opts] = check_inputs('MI',inputs,varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('MI',inputs,varargin{:});
end
possibleOutputs = {'I(A;B)', 'Ish(A;B)',  'Ilin(A;B)', 'coI(A;B)', 'coIsh(A;B)', ...
    'Iss(A)', 'Ic(A;B)', 'Icsh(A;B)', 'Ici(A;B)', 'Icd(A;B)', 'Icdsh(A;B)'};
[isMember, indices] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('MI:invalidOutput', msg);
end


DimsA = size(inputs{1});
DimsB = size(inputs{2});
nTrials = DimsA(end);
if DimsA(end) ~= DimsB(end)
    msg = sprintf('The number of trials for A (%d) and B (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(end),DimsB(end));
    error('MI:InvalidInput', msg);
end

if length(DimsA) > 2 || length(DimsB) > 2
    if length(DimsA) > 2 && length(DimsB) <= 2
        inputs{2} = reshape(inputs{2}, [DimsB(1), 1, DimsB(2)]);
        inputs{2} = repmat(inputs{2}, [1, DimsA(2), 1]);
    elseif length(DimsB) > 2 && length(DimsA) <= 2
        inputs{1} = reshape(inputs{1}, [DimsA(1), 1, DimsA(2)]);
        inputs{1} = repmat(inputs{1}, [1, DimsB(2), 1]);
    end
    opts.timeseries = true;
    nTimepoints = size(inputs{1},2);
else
    nTimepoints = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Prepare Data (binning/reduce dimensions)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~opts.isBinned 
    inputs_b = binning(inputs,opts);
    opts.isBinned = true;
else
    inputs_b = inputs;
end

inputs_1d = inputs_b;
if DimsA(1) > 1
    inputs_1d{1} = reduce_dim(inputs_b{1}, 1);
    if  any(strcmp(reqOutputs,'Hlin(A)')) || any(strcmp(reqOutputs,'Hind(A)')) || any(strcmp(reqOutputs, 'Hind(A|B)'))
        inputs_1d{3} = inputs_b{1};
    end 
end
if DimsB(1) > 1
    inputs_1d{2} = reduce_dim(inputs_b{2}, 1);
end
if DimsA(2:end) ~= DimsB(2:end)
    msg = 'Inconsistent sizes of A and B';
    error('H:inconsistentSizes', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 3: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr = opts.bias;
if opts.computeNulldist == true
    nullDist_opts = opts;
    nullDist_opts.computeNulldist = false;
    MI_nullDist = create_nullDist(inputs_b, reqOutputs, @MI, nullDist_opts);
else
    MI_nullDist = 0;
end

if ~strcmp(corr, 'naive') && ~strcmp(corr, 'bub') && ~strcmp(corr, 'pt')
    [MI_values, MI_naive] = correction(inputs_b, reqOutputs, corr,  @MI, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 4: Compute required Entropies for the requested Outputs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_dependencies = struct( ...
    'I_A_B', {{'H(A)', 'H(A|B)'}}, ...
    'Ish_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Hind(A|B)'}}, ...
    'Ilin_A_B', {{'Hlin(A)', 'H(A|B)', 'Hind(A|B)'}}, ...
    'coI_A_B', {{'H(A)', 'H(A|B)', 'Hlin(A)', 'Hind(A|B)'}}, ...
    'coIsh_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Hlin(A)'}}, ...
    'Iss_A', {{'Hind(A)', 'Hlin(A)'}}, ...
    'Ic_A_B', {{'H(A)', 'H(A|B)', 'Hind(A|B)', 'Hind(A)'}}, ...
    'Icsh_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Hind(A)'}}, ...
    'Ici_A_B', {{'Chi(A)', 'Hind(A)'}}, ...
    'Icd_A_B', {{'H(A)', 'H(A|B)', 'Chi(A)', 'Hind(A|B)'}}, ...
    'Icdsh_A_B', {{'H(A)', 'H(A|B)', 'Hsh(A|B)', 'Chi(A)'}} ...
    );
required_entropies = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'I(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.I_A_B{:}];
        case 'Ish(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ish_A_B{:}];
        case 'Ilin(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ilin_A_B{:}];
        case 'coI(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.coI_A_B{:}];
        case 'coIsh(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.coIsh_A_B{:}];
        case 'Iss(A)'
            required_entropies = [required_entropies, entropy_dependencies.Iss_A{:}];
        case 'Ic(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ic_A_B{:}];
        case 'Icsh(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Icsh_A_B{:}];
        case 'Ici(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Ici_A_B{:}];
        case 'Icd(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Icd_A_B{:}];
        case 'Icdsh(A;B)'
            required_entropies = [required_entropies, entropy_dependencies.Icdsh_A_B{:}];
    end
end
opts_entropies = opts;
opts_entropies.computeNulldist = false;
required_entropies = unique(required_entropies);
if strcmp(corr, 'pt')
    [H_values, H_naive] = H(inputs_b, required_entropies, opts_entropies);
else
    [H_values] = H(inputs_b, required_entropies, opts_entropies);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 5: Compute requested Output Values                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Initialize cell for MI values
MI_values = cell(1, length(reqOutputs));
MI_naive = cell(1, length(reqOutputs));
for t = 1:nTimepoints
    for i = 1:length(indices)
        idx = indices(i);
        switch possibleOutputs{idx}
            case 'I(A;B)'
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                MI_values{i}(1,t) = H_A - H_A_B;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_n - H_A_B_n;
                end
            case 'Ish(A;B)'
                % Ish(A;B) = H(A) - Hind(A|B) + Hsh(A|B) - H(A|B)
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                MI_values{i}(1,t) = H_A - Hind_A_B + Hsh_A_B - H_A_B;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hind_A_B_n = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hsh_A_B_n = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_n - Hind_A_B_n + Hsh_A_B_n - H_A_B_n;
                end
            case 'Ilin(A;B)'
                % I_lin(A;B) = Hlin(A) - Hind(A|B)
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                MI_values{i}(1,t) = Hlin_A - Hind_A_B;
                if strcmp(corr, 'pt')
                    Hlin_A_n = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hind_A_B_n = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = Hlin_A_n - Hind_A_B_n;
                end
            case 'coI(A;B)'
                % coI(A;B) = H(A) - H(A|B) - Hlin(A) + Hind(A|B)
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                MI_values{i}(1,t) = H_A - H_A_B - Hlin_A + Hind_A_B;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hlin_A_n = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_B_n = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_n - H_A_B_n - Hlin_A_n + Hind_A_B_n;
                end
            case 'coIsh(A;B)'
                % coIsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Hlin(A)
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                MI_values{i}(1,t) = H_A + Hsh_A_B - H_A_B - Hlin_A;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hlin_A_n = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hsh_A_B_n = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_n + Hsh_A_B_n - H_A_B_n - Hlin_A_n;
                end
            case 'Iss(A)'
                % Iss(A) = Hind(A) - Hlin(A)
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                MI_values{i}(1,t) = Hind_A - Hlin_A;
                if strcmp(corr, 'pt')
                    Hind_A_n = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    Hlin_A_n = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    MI_naive{i}(1,t) = Hind_A_n - Hlin_A_n;
                end
            case 'Ic(A;B)'% possibleOutputs = { 'H(A|B)', 'Hind(A|B)', 'Hsh(A|B)'};
                % Ic(A;B) = H(A) - H(A|B) + Hind(A|B) - Hind(A)
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                MI_values{i}(1,t) = H_A - H_A_B + Hind_A_B - Hind_A;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_B_n = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hind_A_n = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = H_A_n - H_A_B_n + Hind_A_B_n - Hind_A_n;
                end
            case 'Icsh(A;B)'
                % Icsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Hind(A)
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                MI_values{i}(1,t) = H_A + Hsh_A_B - H_A_B - Hind_A;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_n = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_n = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = H_A_n + Hsh_A_B_n - H_A_B_n - Hind_A_n;
                end
            case 'Ici(A;B)'
                % Ici(A;B) = Chi(A) - Hind(A)
                Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')}(t);
                Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                MI_values{i}(1,t) = Chi_A - Hind_A;
                if strcmp(corr, 'pt')
                    Chi_A_n = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A_n = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = Chi_A_n - Hind_A_n;
                end
            case 'Icd(A;B)'
                % Icd(A;B) = H(A) - H(A|B) - Chi(A) + Hind(A|B)
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')}(t);
                Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                MI_values{i}(1,t) = H_A - H_A_B - Chi_A + Hind_A_B;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A_n = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A_B_n = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_n - H_A_B_n - Chi_A_n + Hind_A_B_n;
                end
            case 'Icdsh(A;B)'
                % Icdsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Chi(A)
                H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')}(t);
                MI_values{i}(1,t) = H_A + Hsh_A_B - H_A_B - Chi_A;
                if strcmp(corr, 'pt')
                    H_A_n = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_n = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_n = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A_n = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    MI_naive{i}(1,t) = H_A_n + Hsh_A_B_n - H_A_B_n - Chi_A_n;
                end
        end
    end
end


