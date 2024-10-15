function [entropies, entropies_naive, entropies_nullDist, prob_dists] = H(inputs, varargin)
%%% *function [entropies, entropies_naive, shuff_all] = H(inputs, outputs, opts)*
%%%
%%% The TE function computes Transfer Entropy (TE) values from time series data.
%%% Transfer entropy from a process A to another process B is the amount of uncertainty
%%% reduced in future values of B by knowing the past values of A given past values of B.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input data A and B. Each cell represents data, where:
%%%             - inputs{1}: First data input (A) with dimensions
%%%                          nDims [X nTimepoints] X nTrials (can be 2 or 3 dimensional)
%%%             - inputs{2}: Second data input (B) with dimensions
%%%                          nDims [X nTimepoints] X nTrials
%%%
%%%             Note: The first dimension (nDims) can accommodate multiple neurons
%%%             concatenated together, allowing for analysis across multiple
%%%             neuronal datasets.
%%%
%%%   - outputs: A cell array of strings specifying which entropies to compute.
%%%              Possible outputs include:
%%%              - 'H(A)'        : Entropy of A.
%%%              - 'H(A|B)'      : Conditional entropy of A given B.
%%%              - 'Hlin(A)'     : Linear entropy of A.
%%%              - 'Hind(A)'     : Induced entropy of A.
%%%              - 'Hind(A|B)'   : Induced conditional entropy of A given B.
%%%              - 'Chi(A)'      : Chi entropy of A.
%%%              - 'Hsh(A)'      : Shuffled entropy of A.
%%%              - 'Hsh(A|B)'    : Shuffled conditional entropy of A given B.
%%%
%%%   - varargin: Optional arguments, passed as a structure. Fields may include:
%%%              - bias:           Specifies the bias correction method to be used.
%%%                                Possible values include:
%%%                                'naive'                      :(default) - No correction applied.
%%%                                'qe', 'le'                   :quadratic/linear extrapolation (need to specify xtrp as number of extrapolations).
%%%                                'ShuffSub'                   :Shuffle Substraction (need to specify shuff as number of shufflings).
%%%                                'qe_ShuffSub', 'le_ShuffSub' :Combination of qe/le and Shuffsub (need to specify shuff and xtrp).
%%%                                'pt'                         :Panzeri-Treves bias correction.
%%%                                'bub'                        :
%%%                                Users can also define their own custom bias correction method
%%%                                (see help for correction.m).
%%%
%%%              - bin_method:     Cell array specifying the binning method.
%%%                                It can have one or two entries:
%%%                                If one entry is provided, it will be applied to both A and B.
%%%                                Possible values include:
%%%                                'none'      : (default) - No binning applied.
%%%                                'eqpop'     : Equal population binning.
%%%                                'eqspace'   : Equal space binning.
%%%                                'threshold' : Binning based on a specified threshold.
%%%                                Users can also define their own custom binning method
%%%                                (see help for binning.m).
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
%%%
%%%
%%% Outputs:
%%%   - entropies: A cell array containing the computed entropy values as specified in the outputs argument.
%%%   - entropies_naive: A cell array containing the naive entropy estimates.
%%%   - shuff_all: A value indicating the results of the shuffling procedure (0 if not performed).
%%%
%%% Note:
%%% Input A and B can represent multiple neurons concatenated along the first dimension.
%%% This means that each neuron can contribute its activity data, allowing the analysis
%%% of interactions between different neurons and their influence on other time series.
%%%
%%% Example:
%%%   To compute the Entropy from two neural populations X1 and X2,
%%%   the function can be called as follows:
%%%   inputs = {X1, X2};
%%%   outputs = {'H(A)', 'H(A|B)'};
%%%   [entropies, entropies_naive, shuff_all] = H(inputs, outputs, opts);


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
        [outputs, opts] = check_inputs('H',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('H',inputs,varargin{:});
end

nVars = length(inputs);
if length(opts.bin_method) < nVars
    opts.bin_method{(end+1):nVars} = opts.bin_method{end};
end
if length(opts.n_bins) < nVars
    opts.n_bins{(end+1):nVars} = opts.n_bins{end};
end

% Check Outputslist
possibleOutputs = {'H(A)', 'H(A|B)', 'Hlin(A)', 'Hind(A)', 'Hind(A|B)', 'Chi(A)','Hsh(A)', 'Hsh(A|B)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('H:invalidOutput', msg);
end

DimsA = size(inputs{1});
DimsB = size(inputs{2});
if DimsA(end) ~= DimsB(end)
    msg = sprintf('The number of trials for A (%d) and B (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(end),DimsB(end));
    error('H:InvalidInput', msg);
end
nTrials = DimsA(end);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Prepare Data (binning/reduce dimensions)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~opts.isbinned 
    inputs_b = binning(inputs,opts);
    opts.isbinned = true;
    for c=1:length(inputs)
        inputs_b{c} = inputs_b{c}+1;
    end
else
    inputs_b = inputs;

end

inputs_1d = inputs_b;
if DimsA(1) > 1
    inputs_1d{1} = reduce_dim(inputs_b{1}, 1);
    if  any(strcmp(outputs,'Hlin(A)')) || any(strcmp(outputs,'Hind(A)')) || any(strcmp(outputs, 'Hind(A|B)'))
        inputs_1d{3} = inputs{1};
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
nullDist_opts.compute_nulldist = false;

entropies_nullDist = 0;
if opts.compute_nulldist
        entropies_nullDist = create_NullDistribution(inputs, outputs, @H, nullDist_opts);
end
if ~strcmp(corr, 'naive')
    [entropies, entropies_naive, entropies_nullDist] = correction(inputs_1d, outputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Step 3.B: Compute required Probability Distributions               %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_distributions = struct( ...
    'H_A', {{'P(A)'}}, ...
    'H_A_B', {{'P(A|B)', 'P(B)'}}, ...
    'Hlin_A', {{'Plin(A)'}}, ...
    'Hind_A', {{'Pind(A)'}}, ...
    'Hind_A_B', {{'Pind(A|B)', 'P(B)'}}, ...
    'Chi_A', {{'P(A)','Pind(A)'}}, ...
    'Hsh_A', {{'Psh(A)'}}, ...
    'Hsh_A_B', {{'Psh(A|B)', 'P(B)'}} ...
    );

required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'H(A)'
            required_distributions = [required_distributions, entropy_distributions.H_A{:}];
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
    end
end
required_distributions = unique(required_distributions);
prob_dists = prob_estimator(inputs_1d, required_distributions, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Step 4.B: Compute requested Entropies                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropies = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'H(A)'
            P_A = prob_dists{strcmp(required_distributions, 'P(A)')};
            P_lin_log = P_A .* log2(P_A);
            P_lin_log(isnan(P_lin_log)) = 0;
            entropies{i} = -sum(P_lin_log(:));
        case 'H(A|B)'
            P_AB = prob_dists{strcmp(required_distributions, 'P(A|B)')};
            P_B = prob_dists{strcmp(required_distributions, 'P(B)')};
            P_Bext = repmat(P_B,1, size(P_AB,1));
            P_lin_log = P_Bext' .* P_AB .* log2(P_AB);
            P_lin_log(isnan(P_lin_log)) = 0;
            entropies{i} = -sum(P_lin_log(:));
        case 'Hlin(A)'
            P_lin = prob_dists{strcmp(required_distributions, 'Plin(A)')};
            P_lin_log = P_lin .* log2(P_lin);
            P_lin_log(isnan(P_lin_log)) = 0;  % Handle 0*log(0) cases
            entropies{i} = -sum(P_lin_log(:));
        case 'Hind(A)'
            P_indA = prob_dists{strcmp(required_distributions, 'Pind(A)')};
            P_lin_log = P_indA .* log2(P_indA);
            P_lin_log(isnan(P_lin_log)) = 0;
            entropies{i} = -sum(P_lin_log(:));
        case 'Hind(A|B)'
            P_indAB = prob_dists{strcmp(required_distributions, 'Pind(A|B)')};
            P_B = prob_dists{strcmp(required_distributions, 'P(B)')};
            P_lin_log = P_B' .* P_indAB .* log2(P_indAB);
            P_lin_log(isnan(P_lin_log)) = 0;
            entropies{i} = -sum(P_lin_log(:));
        case 'Chi(A)'
            P_A = prob_dists{strcmp(required_distributions, 'P(A)')};
            P_indA = prob_dists{strcmp(required_distributions, 'Pind(A)')};
            if size(P_A,1)<size(P_indA,1)
               P_A = [P_A; zeros(size(P_indA,1)-size(P_A,1))];
            end
            P_lin_log = P_A .* log2(P_indA);
            P_lin_log(isnan(P_lin_log)) = 0;
            P_lin_log(isinf(P_lin_log)) = 0;
            entropies{i} = -sum(P_lin_log(:));
        case 'Hsh(A)'
            P_shA = prob_dists{strcmp(required_distributions, 'Psh(A)')};
            P_lin_log = P_shA .* log2(P_shA);
            P_lin_log(isnan(P_lin_log)) = 0;
            entropies{i} = -sum(P_lin_log(:));
        case 'Hsh(A|B)'
            P_shAB = prob_dists{strcmp(required_distributions, 'Psh(A|B)')};
            P_B = prob_dists{strcmp(required_distributions, 'P(B)')};
            P_lin_log = P_B' .* P_shAB .* log2(P_shAB);
            P_lin_log(isnan(P_lin_log)) = 0;
            entropies{i} = -sum(P_lin_log(:));
    end
end
entropies_naive = entropies;
end

