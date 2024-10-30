function [cMI_values, cMI_naive, cMI_nullDist] =  cMI(inputs, varargin)
% *function [cMI_values, cMI_naive, cMI_nullDist] =  cMI(inputs, reqOutputs, opts)*
%
% This function calculates conditional mutual information (cMI) 
%
% Inputs:
%   - inputs: A cell array containing the data:
%             - inputs{1}: First input data (A) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             - inputs{2}: Second input data (B) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             - inputs{3}: Condition input data (C) with dimensions
%                          nDims X (nTimepoints X) nTrials
%             -> In cases where the input is provided as a time series, the reqOutputs 
%                will be computed for each time point, resulting in outputs that are 
%                also represented as time series
%
%   - reqOutputs: A cell array of strings specifying which entropies or cMI measures to compute.
%               - 'I(A;B|C)'    : Conditional Mutual Information I(A;B|C)
%               - 'Ish(A;B|C)'  : Shuffled Conditional MI Ish(A;B|C)
%               - 'Ilin(A;B|C)' : Linear Conditional MI Ilin(A;B|C)
%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - bias:             Specifies the bias correction method to be used.
%                                  'naive'                      :(default) - No correction applied.
%                                  'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%                                  'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%                                  'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%                                  'pt'                         :Panzeri-Treves bias correction (Panzeri and Treves 1996).
%                                  'bub'                        :best upper bound(Paninsky, 2003)
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
%   - cMI_values: A cell array containing the computed cMI values as specified in the reqOutputs argument.
%   - cMI_naive: A cell array containing the naive cMI estimates.
%   - cMI_shuff_all: Results of the null distribution computation (0 if not performed).
%
% Example:
% Suppose we have two time series of neural activity X1 and X2 and a
% time series Y. To compute the conditional mutual information I(X1;X2|Y)
% over time we have to call cMI as follows:
%   cMI_values = cMI({X1, X2, Y}, {'I(A;B|C)'}, opts);
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
        [inputs, reqOutputs, opts] = check_inputs('cMI',inputs,varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('cMI',inputs,varargin{:});
end

possibleOutputs = {'I(A;B|C)', 'Ish(A;B|C)', 'Ilin(A;B|C)'};
[isMember, indices] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('cMI:invalidOutput', msg);
end
DimsA = size(inputs{1});
if length(DimsA) > 2
    nTimepoints = DimsA(2);
else
    nTimepoints = 1;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%        Step 2: Compute required Entropies for the requested Outputs           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_dependencies = struct( ...
    'I_AB_C', {{'H(A|C)', 'H(A|B,C)'}} ...
    );

required_entropies = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'I(A;B|C)'
            required_entropies = [required_entropies, entropy_dependencies.I_AB_C{:}];
            % case 'Ish(A;B|C)'
            %     required_entropies = [required_entropies, entropy_dependencies.Ish_AB_C{:}];
            % case 'Ilin(A;B|C)'
            %     required_entropies = [required_entropies, entropy_dependencies.Ilin_AB_C{:}];
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
        case 'H(A|C)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs{1}, inputs{3}}, {'H(A|B)'}, opts_entropy);          
        case 'H(A|B,C)'
            [H_values{i}, H_naive{i}, H_shuff_all{i}] = H({inputs{1}, cat(1,inputs{2}, inputs{3})}, {'H(A|B)'}, opts_entropy);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%H_naiv
%                  Step 3: Compute requested Output Values                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Initialize cell for MI values
cMI_values = cell(nTimepoints, length(reqOutputs));
cMI_naive = cell(nTimepoints, length(reqOutputs));
cMI_shuff_all = cell(nTimepoints, length(reqOutputs));
for t = 1:nTimepoints    
    if iscell(H_shuff_all)
        H_t_shuff_all = H_shuff_all(t, :);
    end
    for i = 1:length(indices)
        idx = indices(i);
        nOut = nargout;
        switch possibleOutputs{idx}
            case 'I(A;B|C)'
                % I(A;B|C) = H(A|C) - H(A|B,C)
                %____________________________________________________________________________________________%
                if strcmp(opts.bias, 'shuffSub')
                    H_AC_naive = H_naive{strcmp(required_entropies, 'H(A|C)')}(t);
                    H_ABC_naive = H_naive{strcmp(required_entropies, 'H(A|B,C)')}(t);
                    cMI_naive{t, i} =  H_AC_naive - H_ABC_naive;
                    H_AC_shuff_all = cell2mat(H_t_shuff_all{strcmp(required_entropies, 'H(A|C)')});
                    H_ABC_shuff_all = cell2mat(H_t_shuff_all{strcmp(required_entropies, 'H(A|B,C)')});
                    for shuffIdx = 1:opts.shuff
                        H_AC_shuff = H_AC_shuff_all(shuffIdx);
                        H_ABC_shuff = H_ABC_shuff_all(shuffIdx);
                        cMI_shuff_all{t, i} =  H_AC_shuff - H_ABC_shuff;
                    end
                    cMI_values{t, i} = cMI_naive{t, i} - mean(cMI_shuff_all{t, i});
                    nOut = 1;
                else
                    H_AC = H_values{strcmp(required_entropies, 'H(A|C)')}(t);
                    H_ABC = H_values{strcmp(required_entropies, 'H(A|B,C)')}(t);
                    cMI_values{t, i} = H_AC -  H_ABC;
                end
                %____________________________________________________________________________________________%
                if nOut > 1
                    H_AC_naive = H_naive{strcmp(required_entropies, 'H(A|C)')}(t);
                    H_ABC_naive = H_naive{strcmp(required_entropies, 'H(A|B,C)')}(t);
                    cMI_naive{t, i} =  H_AC_naive - H_ABC_naive;
                    H_AC_shuff_all = cell2mat(H_t_shuff_all{strcmp(required_entropies, 'H(A|C)')});
                    H_ABC_shuff_all = cell2mat(H_t_shuff_all{strcmp(required_entropies, 'H(A|B,C)')});
                end
        end
    end
end
end

