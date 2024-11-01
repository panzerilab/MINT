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
if length(DimsA) > 2
    nTimepointsA = DimsA(2);
else
    nTimepointsA = 1;
end
DimsB = size(inputs{2});
if length(DimsB) > 2
    nTimepointsB = DimsB(2);
else
    nTimepointsB = 1;
end
nTimepoints = max(nTimepointsA, nTimepointsB);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 2: Compute required Entropies for the requested Outputs           %
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
required_entropies = unique(required_entropies);
corr = opts.bias;
opts_entropies = opts;
opts_entropies.computeNulldist = false;
%opts_entropies.isBinned = false;
if strcmp(opts.bias, 'shuffCorr')
    [MI_values, MI_naive, MI_nullDist] = correction(inputs, reqOutputs, corr, @MI, opts_entropies);
    return
else 
    MI_nullDist = 0;
end 
[H_values, H_naive, H_shuff_all] = H(inputs, required_entropies, opts_entropies);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 3: Compute requested Output Values                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize cell for MI values
MI_values = cell(1, length(reqOutputs));
MI_naive = cell(1, length(reqOutputs));

for t = 1:nTimepoints
    if iscell(H_shuff_all)
        for i = 1:length(H_values)
             H_t_shuff_all{i} = H_shuff_all{i}(:, t);
        end
    end
    for i = 1:length(indices)
        idx = indices(i);
        nOut = nargout;
        switch possibleOutputs{idx}
            case 'I(A;B)'
                % I(A;B) = H(A) - H(A|B)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive;

                    % Loop through shuffled values for this timepoint
                    for shuffIdx = 1:opts.shuff
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = H_A_naive - H_A_B_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_values{i}(1,t) = H_A - H_A_B;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive;

                end
            case 'Ish(A;B)'
                % Ish(A;B) = H(A) - Hind(A|B) + Hsh(A|B) - H(A|B)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - Hind_A_B_naive + Hsh_A_B_naive - H_A_B_naive;
                    for shuffIdx = 1:opts.shuff
                        Hind_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        Hsh_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = H_A_naive - Hind_A_B_shuff + Hsh_A_B_shuff - H_A_B_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_values{i}(1,t) = H_A - Hind_A_B + Hsh_A_B - H_A_B;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - Hind_A_B_naive + Hsh_A_B_naive - H_A_B_naive;

                end

            case 'Ilin(A;B)'
                % I_lin(A;B) = Hlin(A) - Hind(A|B)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = Hlin_A_naive - Hind_A_B_naive;
                    for shuffIdx = 1:opts.shuff
                        Hind_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = Hlin_A_naive - Hind_A_B_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_values{i}(1,t) = Hlin_A - Hind_A_B;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = Hlin_A_naive - Hind_A_B_naive;

                end

            case 'coI(A;B)'
                % coI(A;B) = H(A) - H(A|B) - Hlin(A) + Hind(A|B)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive - Hlin_A_naive + Hind_A_B_naive;
                    for shuffIdx = 1:opts.shuff
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Hind_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = H_A_naive - H_A_B_shuff - Hlin_A_naive + Hind_A_B_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_values{i}(1,t) = H_A - H_A_B - Hlin_A + Hind_A_B;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive - Hlin_A_naive + Hind_A_B_naive;

                end
            case 'coIsh(A;B)'
                % coIsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Hlin(A)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hlin_A_naive;
                    for shuffIdx = 1:opts.shuff
                        Hsh_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Hlin_A_naive;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                    Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    MI_values{i}(1,t) = H_A + Hsh_A_B - H_A_B - Hlin_A;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hlin_A_naive;

                end
            case 'Iss(A)'
                % Iss(A) = Hind(A) - Hlin(A)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    MI_naive{i}(1,t) = Hind_A_naive - Hlin_A_naive;
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff= H_t_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = Hind_A_shuff - Hlin_A_naive;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                    Hlin_A = H_values{strcmp(required_entropies, 'Hlin(A)')}(t);
                    MI_values{i}(1,t) = Hind_A - Hlin_A;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    Hlin_A_naive = H_naive{strcmp(required_entropies, 'Hlin(A)')}(t);
                    MI_naive{i}(1,t) = Hind_A_naive - Hlin_A_naive;
                end
            case 'Ic(A;B)'% possibleOutputs = { 'H(A|B)', 'Hind(A|B)', 'Hsh(A|B)'};
                % Ic(A;B) = H(A) - H(A|B) + Hind(A|B) - Hind(A)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive + Hind_A_B_naive - Hind_A_naive;
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Hind_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = H_A_naive - H_A_B_shuff + Hind_A_B_shuff - Hind_A_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_values{i}(1,t) = H_A - H_A_B + Hind_A_B - Hind_A;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive + Hind_A_B_naive - Hind_A_naive;
                end
            case 'Icsh(A;B)'
                % Icsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Hind(A)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hind_A_naive;
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        Hsh_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Hind_A_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_values{i}(1,t) = H_A + Hsh_A_B - H_A_B - Hind_A;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Hind_A_naive;
                end
            case 'Ici(A;B)'
                % Ici(A;B) = Chi(A) - Hind(A)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = Chi_A_naive - Hind_A_naive;
                    for shuffIdx = 1:opts.shuff
                        Hind_A_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A)')}(shuffIdx);
                        Chi_A_shuff = H_t_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = Chi_A_shuff - Hind_A_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A = H_values{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_values{i}(1,t) = Chi_A - Hind_A;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A_naive = H_naive{strcmp(required_entropies, 'Hind(A)')}(t);
                    MI_naive{i}(1,t) = Chi_A_naive - Hind_A_naive;
                end
            case 'Icd(A;B)'
                % Icd(A;B) = H(A) - H(A|B) - Chi(A) + Hind(A|B)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive - Chi_A_naive + Hind_A_B_naive;
                    for shuffIdx = 1:opts.shuff
                        Chi_A_shuff = H_t_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Hind_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hind(A|B)')}(shuffIdx);
                        MI_shuff_all(shuffIdx)= H_A_naive - H_A_B_shuff - Chi_A_shuff+ Hind_A_B_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A_B = H_values{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_values{i}(1,t) = H_A - H_A_B - Chi_A + Hind_A_B;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    Hind_A_B_naive = H_naive{strcmp(required_entropies, 'Hind(A|B)')}(t);
                    MI_naive{i}(1,t) = H_A_naive - H_A_B_naive - Chi_A_naive + Hind_A_B_naive;
                end
            case 'Icdsh(A;B)'
                % Icdsh(A;B) = H(A) + Hsh(A|B) - H(A|B) - Chi(A)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Chi_A_naive;
                    for shuffIdx = 1:opts.shuff
                        Hsh_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'Hsh(A|B)')}(shuffIdx);
                        H_A_B_shuff = H_t_shuff_all{strcmp(required_entropies, 'H(A|B)')}(shuffIdx);
                        Chi_A_shuff = H_t_shuff_all{strcmp(required_entropies, 'Chi(A)')}(shuffIdx);
                        MI_shuff_all(shuffIdx) = H_A_naive + Hsh_A_B_shuff - H_A_B_shuff - Chi_A_shuff;
                    end
                    MI_values{i}(1,t) = MI_naive{i}(1,t) - mean(MI_shuff_all);
                    nOut = 1;
                else
                    H_A = H_values{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B = H_values{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B = H_values{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A = H_values{strcmp(required_entropies, 'Chi(A)')}(t);
                    MI_values{i}(1,t) = H_A + Hsh_A_B - H_A_B - Chi_A;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_A_naive = H_naive{strcmp(required_entropies, 'H(A)')}(t);
                    Hsh_A_B_naive = H_naive{strcmp(required_entropies, 'Hsh(A|B)')}(t);
                    H_A_B_naive = H_naive{strcmp(required_entropies, 'H(A|B)')}(t);
                    Chi_A_naive = H_naive{strcmp(required_entropies, 'Chi(A)')}(t);
                    MI_naive{i}(1,t) = H_A_naive + Hsh_A_B_naive - H_A_B_naive - Chi_A_naive;
                end
        end

        if opts.computeNulldist
            nullDist_opts = opts;
            nullDist_opts.computeNulldist = false;
            MI_nullDist = create_nullDist(inputs, reqOutputs, @MI, nullDist_opts);
        end
    end
end


