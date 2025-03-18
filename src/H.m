function [entropies, entropies_plugin, entropies_nullDist, prob_dists] = H(inputs, varargin)
% *function [entropies, entropies_plugin, entropies_nullDist, prob_dists] = H(inputs, reqOutputs, opts)*
% H - Calculate Entropy (H) and related information-theoretic quantities
%
% This function calculates entropy (H) and other related measures based on the 
% provided inputs, reqOutputs, and optional parameters.
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
%               - 'H(A)'        : Entropy of A.
%               - 'H(B)'        : Entropy of B.
%               - 'H(A|B)'      : Conditional entropy of A given B.
%               - 'Hlin(A)'     : Linear entropy of A.
%               - 'Hind(A)'     : Independent entropy of A.
%               - 'Hind(A|B)'   : Independent conditional entropy of A given B.
%               - 'Chi(A)'      : Chi entropy of A.
%               - 'Hsh(A)'      : Shuffled entropy of A.
%               - 'Hsh(A|B)'    : Shuffled conditional entropy of A given B.
%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - bias:             Specifies the bias correction method to be used.
%                                  'plugin'                     :(default) - No correction applied.
%                                  'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%                                  'nsb'                        :correction using the NSB algorithm (Nemenman, Bialek and van Steveninck, 2019) 
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
%                                  'setToZero'   : Sets NaN values to zero.
%                                  'error'       : (default) Throws an error if NaN values are detected.
%
% Outputs:
%   - entropies: A cell array containing the computed entropy values as specified in the reqOutputs argument.
%   - entropies_plugin: A cell array containing the plugin entropy estimates.
%   - entropies_nullDist: Results of the null distribution computation (0 if not performed).
%   - prob_dists: A cell array containing the estimated probability distributions used in entropy calculations.
%
% EXAMPLE
% Suppose we have two time series of groups of neurons X1 and X2
% and two time series of groups of neurons Y1 and Y2.
% (Structure of X1, X2, Y1, Y2 is nNeurons x nTrials)
%
% We can structure our inputs as follows:
% X = cat(1, X1, X2);  % Concatenates X1 and X2 along the first dimension (neurons)
% Y = cat(1, Y1, Y2);  % Concatenates Y1 and Y2 along the first dimension (neurons)
%
% To compute the conditioned Entropy and the Shuffled Entropy of X (X|Y), the function can be called as:
% [entropies, entropies_plugin] = H({X, Y}, {'H(A)', 'H(A|B)'}, opts);
%
% Here, 'opts' represents additional options you may want to include (see varargin options).

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
        [inputs, reqOutputs, opts] = check_inputs('H',inputs,varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('H',inputs,varargin{:});
end

nVars = length(inputs);
if length(opts.bin_method) < nVars
    opts.bin_method{(end+1):nVars} = opts.bin_method{end};
end
if length(opts.n_bins) < nVars
    opts.n_bins{(end+1):nVars} = opts.n_bins{end};
end

% Check Outputslist
possibleOutputs = {'H(A)', 'H(B)','H(A|B)', 'Hlin(A)', 'Hind(A)', 'Hind(A|B)', 'Chi(A)','Hsh(A)', 'Hsh(A|B)', 'Hnsb(A)', 'Hnsb(A,B)', 'Hnsb(B)'};
[isMember, indices] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('H:invalidOutput', msg);
end

DimsA = size(inputs{1});
DimsB = size(inputs{2});
nTrials = DimsA(end);
if DimsA(end) ~= DimsB(end)
    msg = sprintf('The number of trials for A (%d) and B (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(end),DimsB(end));
    error('H:InvalidInput', msg);
end

if length(DimsA) > 2 || length(DimsB) > 2
    if length(DimsA) > 2 && length(DimsB) <= 2
        inputs{2} = reshape(inputs{2}, [DimsB(1), 1, DimsB(2)]);
        inputs{2} = repmat(inputs{2}, [1, DimsA(2), 1]);
    elseif length(DimsB) > 2 && length(DimsA) <= 2
        inputs{1} = reshape(inputs{1}, [DimsA(1), 1, DimsA(2)]);
        inputs{1} = repmat(inputs{1}, [1, DimsB(2), 1]);
    end
    nTimepoints = (DimsA(2));
    opts.timeseries = true;
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
%                    Step 3.A: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr = opts.bias;
corefunc = @H;
nullDist_opts = opts;
nullDist_opts.computeNulldist = false;

if opts.computeNulldist == true
        entropies_nullDist = create_nullDist(inputs_b, reqOutputs, @H, nullDist_opts);
else 
    entropies_nullDist = 0;
end 

if ~strcmp(corr, 'plugin') &&  ~strcmp(corr, 'bub') &&  ~strcmp(corr, 'pt')
    [entropies, entropies_plugin, entropies_shuffAll] = correction(inputs_1d, reqOutputs, corr, corefunc, opts);
    if ~iscell(entropies_nullDist)
        entropies_nullDist = entropies_shuffAll;
    end 
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Step 3.B: Compute required Probability Distributions               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_distributions = struct( ...
    'H_A', {{'P(A)'}}, ...
    'H_B', {{'P(B)'}}, ...
    'H_A_B', {{'P(A|B)', 'P(B)', 'P(A)'}}, ...
    'Hlin_A', {{'Plin(A)'}}, ...
    'Hind_A', {{'Pind(A)'}}, ...
    'Hind_A_B', {{'Pind(A|B)', 'P(B)'}}, ...
    'Chi_A', {{'P(A)','Pind(A)'}}, ...
    'Hsh_A', {{'Psh(A)'}}, ...
    'Hsh_A_B', {{'Psh(A)','Psh(A|B)', 'P(B)'}}, ...
    'Hnsb_A', {{'P(A)'}}, ...
    'Hnsb_B', {{'P(B)'}}, ...
    'Hnsb_AB', {{'P(A,B)'}} ...
    );

required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'H(A)'
            required_distributions = [required_distributions, entropy_distributions.H_A{:}];
        case 'H(B)'
            required_distributions = [required_distributions, entropy_distributions.H_B{:}];
        case 'H(A|B)'
            required_distributions = [required_distributions, entropy_distributions.H_A_B{:}];
        case 'Hlin(A)'
            required_distributions = [required_distributions, entropy_distributions.Hlin_A{:}];
        case 'Hind(A)'
            required_distributions = [required_distributions, entropy_distributions.Hind_A{:}];
        case 'Hind(A|B)'
            required_distributions = [required_distributions, entropy_distributions.Hind_A_B{:}];
        case 'Chi(A)'
            required_distributions = [required_distributions, entropy_distributions.Chi_A{:}];
        case 'Hsh(A)'
            required_distributions = [required_distributions, entropy_distributions.Hsh_A{:}];
        case 'Hsh(A|B)'
            required_distributions = [required_distributions, entropy_distributions.Hsh_A_B{:}];
        case 'Hnsb(A)'
            required_distributions = [required_distributions, entropy_distributions.Hnsb_A{:}];
        case 'Hnsb(B)'
            required_distributions = [required_distributions, entropy_distributions.Hnsb_B{:}];
        case 'Hnsb(A,B)'
            required_distributions = [required_distributions, entropy_distributions.Hnsb_AB{:}];
    end
end
required_distributions = unique(required_distributions);
prob_dists = prob_estimator(inputs_1d, required_distributions, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Step 4.B: Compute requested Entropies                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropies = cell(1, length(reqOutputs));
entropies_plugin = cell(1, length(reqOutputs));

for t = 1:nTimepoints
    for i = 1:length(indices)
        idx = indices(i);
        switch possibleOutputs{idx}
            case 'H(A)'
                P_A = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                P_lin_log = P_A .* log2(P_A);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                if strcmp(opts.bias, 'bub')
                    bias = bub(nTrials * P_A);                           
                elseif strcmp(opts.bias, 'pt')
                    bias = pt(inputs_1d{1}, length(unique(inputs_1d{1})), nTrials);
                else 
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;   
            case 'H(B)'
                P_B = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                P_lin_log = P_B .* log2(P_B);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                if strcmp(opts.bias, 'bub')
                    bias = bub(nTrials * P_B);
                elseif strcmp(opts.bias, 'pt')
                    bias = pt(inputs_1d{2}, length(unique(inputs_1d{2})), nTrials);
                else 
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;   
            case 'H(A|B)'
                P_AB = prob_dists{t, strcmp(required_distributions, 'P(A|B)')};
                P_B = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                P_Bext = repmat(P_B, 1, size(P_AB, 1));
                P_lin_log = P_Bext' .* P_AB .* log2(P_AB);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                if strcmp(opts.bias, 'bub')
                    P_A = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                    pAB = P_AB .* P_B';
                    biasA  = bub(nTrials*P_A(:));
                    biasAB = bub(nTrials*pAB(:));
                    bias = biasAB - biasA;
                elseif strcmp(opts.bias, 'pt')
                    bias = 0;
                    uniqueB = unique(inputs_1d{2});
                    nb = length(uniqueB);
                    for b_i=1:nb
                        A_tmp = inputs_1d{1}(inputs_1d{2} == uniqueB(b_i));
                        bias = bias + pt(A_tmp, length(unique(A_tmp)), nTrials);
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
            case 'Hlin(A)'
                P_lin = prob_dists{t, strcmp(required_distributions, 'Plin(A)')};
                P_lin_log = P_lin .* log2(P_lin);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                if strcmp(opts.bias, 'bub')
                    bias = 0;
                    for row=1:size(plin,1)
                        pl = P_lin(row,:);
                        bias =  bias + bub(nTrials*pl);
                    end
                elseif strcmp(opts.bias, 'pt')
                    bias = 0;
                    for row=1:size(inputs{1},1)
                        Arow = inputs{1}(row,:);
                        nArow = length(unique(Arow));
                        bias =  bias + pt(Arow, nArow, nTrials);
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
            case 'Hind(A)'
                P_indA = prob_dists{t, strcmp(required_distributions, 'Pind(A)')};
                P_lin_log = P_indA .* log2(P_indA);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                % bias = 0;
                if strcmp(opts.bias, 'pt')
                    bias = 0;
                    uniqueB = unique(inputs_1d{2});
                    nb = length(uniqueB);
                    for b_i=1:nb
                        A_tmp = inputs_1d{1}(inputs_1d{2} == uniqueB(b_i));
                        bias = bias + pt(A_tmp, length(unique(A_tmp)), nTrials);
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
            case 'Hind(A|B)'
                P_indAB = prob_dists{t, strcmp(required_distributions, 'Pind(A|B)')};
                P_B = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                P_lin_log = P_B' .* P_indAB .* log2(P_indAB);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                bias = 0;
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
            case 'Chi(A)'
                P_A = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                P_indA = prob_dists{t, strcmp(required_distributions, 'Pind(A)')};
                if size(P_A, 1) < size(P_indA, 1)
                   P_A = [P_A; zeros(size(P_indA, 1) - size(P_A, 1), 1)];
                elseif size(P_A, 1) > size(P_indA, 1)
                    P_indA = [P_indA; zeros(size(P_A, 1) - size(P_indA, 1), 1)];
                end
                P_lin_log = P_A .* log2(P_indA);
                P_lin_log(isnan(P_lin_log)) = 0;
                P_lin_log(isinf(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                bias = 0;
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
            case 'Hsh(A)'
                P_shA = prob_dists{t, strcmp(required_distributions, 'Psh(A)')};
                P_lin_log = P_shA .* log2(P_shA);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                if strcmp(opts.bias, 'bub')
                    bias = bub(nTrials*P_shA);
                    entropies_plugin{i}(1,t) = entropies{i}(1,t);
                    entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
                elseif strcmp(opts.bias, 'pt')
                    bias = pt(inputs_1d{1}, length(unique(inputs_1d{1})), nTrials);
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
            case 'Hsh(A|B)'
                P_shAB = prob_dists{t, strcmp(required_distributions, 'Psh(A|B)')};
                P_B = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                P_lin_log = P_B' .* P_shAB .* log2(P_shAB);
                P_lin_log(isnan(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = -sum(P_lin_log(:));
                if strcmp(opts.bias, 'bub')
                    P_shA = prob_dists{t, strcmp(required_distributions, 'Psh(A)')};
                    PshAB = P_shAB .* P_B;
                    biasA  = bub(nTrials*P_shA);
                    biasAB = bub(nTrials*PshAB);
                    bias = biasAB - biasA;
                    entropies_plugin{i}(1,t) = entropies{i}(1,t);
                    entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;
                elseif strcmp(opts.bias, 'pt')
                    bias = 0;
                    uniqueB = unique(inputs_1d{2});
                    nb = length(uniqueB);
                    for b_i=1:nb
                        A_tmp = inputs_1d{1}(inputs_1d{2} == uniqueB(b_i));
                        bias = bias + pt(A_tmp, length(unique(A_tmp)), nTrials);
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;         
            case 'Hnsb(A)'
                qfun = 1;
                precision = .1;
                P_A = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                nb  = P_A' * nTrials;
                K   = length(nb);
                nxa = nb(nb>0);
                kxb = ones(size(nxa));
                
                [Sb_nsb, ~, ~, ~, ~, S_ml,~]=find_nsb_entropy (kxb, nxa, K, precision,qfun);
                entropies_plugin{i}(1,t) = S_ml;
                entropies{i}(1,t) = Sb_nsb;
            case 'Hnsb(B)'
                qfun = 1;
                precision = .1;
                P_B = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                nb  = P_B'*nTrials;
                K   = length(nb);
                nxb = nb(nb>0);
                kxb = ones(size(nxb));
                
                [Sb_nsb, ~, ~, ~, ~, S_ml,~]=find_nsb_entropy (kxb, nxb, K, precision,qfun);
                entropies_plugin{i}(1,t) = S_ml;
                entropies{i}(1,t) = Sb_nsb;
            case 'Hnsb(A,B)'
                qfun = 1;
                precision = .1;
                P_AB = prob_dists{t, strcmp(required_distributions, 'P(A,B)')};
                nab  = P_AB(:)'*nTrials;
                K = length(nab);
                nxab = nab(nab>0);
                kxab = ones(size(nxab));
 
                [Sab_nsb, ~, ~, ~, ~, S_ml,~]=find_nsb_entropy (kxab, nxab, K, precision,qfun);
                entropies_plugin{i}(1,t) = S_ml;
                entropies{i}(1,t) = Sab_nsb;

        end
    end
end
end
%
