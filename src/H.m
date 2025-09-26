function [entropies, entropies_plugin, entropies_nullDist, prob_dists] = H(inputs, varargin)
% *function [entropies, entropies_plugin, entropies_nullDist, prob_dists] = H(inputs, reqOutputs, opts)*
% H - Calculate Entropy (H) and related information-theoretic quantities
%
% Works with full-support probabilities from prob_estimator:
%   P(A), P(B)         -> column vectors
%   P(A,B), P(A|B)     -> |A| x |B|
%   Pind(A), Pind(A|B) -> aligned to |A|
%   Plin(A)            -> per-dimension marginals of A (rows=dims, cols=bins per dim; zero-padded)
%   Psh(A), Psh(A|B)   -> aligned to P(A) / P(A|B)
%
% Copyright (C) 2024
% GPLv3+

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check Inputs, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked') && opts.isChecked
        reqOutputs = varargin{1};
    else
        [inputs, reqOutputs, opts] = check_inputs('H', inputs, varargin{:});
    end
else
    [inputs, reqOutputs, opts] = check_inputs('H', inputs, varargin{:});
end

nVars = length(inputs);
if length(opts.bin_method) < nVars
    opts.bin_method(end+1:nVars) = {opts.bin_method{end}};
end
if length(opts.n_bins) < nVars
    opts.n_bins(end+1:nVars) = {opts.n_bins{end}};
end

% Allowed outputs
possibleOutputs = {'H(A)','H(B)','H(A|B)','Hlin(A)','Hind(A)','Hind(A|B)', ...
                   'Chi(A)','Hsh(A)','Hsh(A|B)','Hnsb(A)','Hnsb(A,B)','Hnsb(B)'};
[isMember, indices] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    error('H:invalidOutput','Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
end

DimsA = size(inputs{1});
DimsB = size(inputs{2});
nTrials = DimsA(end);
if DimsA(end) ~= DimsB(end)
    error('H:InvalidInput', ...
          'The number of trials for A (%d) and B (%d) are not consistent.', DimsA(end), DimsB(end));
end

% Time-series harmonization
if length(DimsA) > 2 || length(DimsB) > 2
    if length(DimsA) > 2 && length(DimsB) <= 2
        inputs{2} = reshape(inputs{2}, [DimsB(1), 1, DimsB(2)]);
        inputs{2} = repmat(inputs{2}, [1, DimsA(2), 1]);
    elseif length(DimsB) > 2 && length(DimsA) <= 2
        inputs{1} = reshape(inputs{1}, [DimsA(1), 1, DimsA(2)]);
        inputs{1} = repmat(inputs{1}, [1, DimsB(2), 1]);
    end
    nTimepoints = DimsA(2);
    opts.timeseries = true;
else
    nTimepoints = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Prepare Data (binning/reduce dimensions)                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if ~opts.isBinned
    inputs_b = binning(inputs, opts);
    opts.isBinned = true;
else
    inputs_b = inputs;
end

inputs_1d = inputs_b;
if DimsA(1) > 1
    inputs_1d{1} = reduce_dim(inputs_b{1}, 1);             % 1D index of A patterns
    % Keep full A for per-dim marginals / independent models
    if any(strcmp(reqOutputs,'Hlin(A)')) || any(strcmp(reqOutputs,'Hind(A)')) || any(strcmp(reqOutputs,'Hind(A|B)'))
        inputs_1d{3} = inputs_b{1};                         % stash full A in slot 3 (as used by prob_estimator)
    end
end
if DimsB(1) > 1
    inputs_1d{2} = reduce_dim(inputs_b{2}, 1);
end

% Safe size check for trailing dims
if (numel(DimsA) > 1 || numel(DimsB) > 1) && ~isequal(DimsA(2:end), DimsB(2:end))
    error('H:inconsistentSizes','Inconsistent sizes of A and B');
end

% Effective bin counts (used by bias corrections)
if ~strcmp(opts.bin_method{1}, 'none')
    nbinsA = DimsA(1) * opts.n_bins{1};
else
    nbinsA = numel(unique(inputs_1d{1}(:)));
end
if ~strcmp(opts.bin_method{2}, 'none')
    nbinsB = DimsB(1) * opts.n_bins{2};
else
    nbinsB = numel(unique(inputs_1d{2}(:)));
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 3.A: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr = opts.bias;
corefunc = @H;
nullDist_opts = opts;
nullDist_opts.computeNulldist = false;

if opts.computeNulldist
    entropies_nullDist = create_nullDist(inputs_b, reqOutputs, @H, nullDist_opts);
else
    entropies_nullDist = 0;
end

% Delegate to correction() for non-plugin/PT/BUB
if ~strcmp(corr, 'plugin') && ~strcmp(corr, 'bub') && ~strcmp(corr, 'pt')
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
    'H_A',       {{'P(A)'}}, ...
    'H_B',       {{'P(B)'}}, ...
    'H_A_B',     {{'P(A|B)', 'P(B)', 'P(A)'}}, ...
    'Hlin_A',    {{'Plin(A)'}}, ...
    'Hind_A',    {{'Pind(A)'}}, ...
    'Hind_A_B',  {{'Pind(A|B)', 'P(B)'}}, ...
    'Chi_A',     {{'P(A)', 'Pind(A)'}}, ...
    'Hsh_A',     {{'Psh(A)'}}, ...
    'Hsh_A_B',   {{'Psh(A)', 'Psh(A|B)', 'P(B)'}}, ...
    'Hnsb_A',    {{'P(A)'}}, ...
    'Hnsb_B',    {{'P(B)'}}, ...
    'Hnsb_AB',   {{'P(A,B)'}} );

required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'H(A)'        , required_distributions = [required_distributions, entropy_distributions.H_A{:}];
        case 'H(B)'        , required_distributions = [required_distributions, entropy_distributions.H_B{:}];
        case 'H(A|B)'      , required_distributions = [required_distributions, entropy_distributions.H_A_B{:}];
        case 'Hlin(A)'     , required_distributions = [required_distributions, entropy_distributions.Hlin_A{:}];
        case 'Hind(A)'     , required_distributions = [required_distributions, entropy_distributions.Hind_A{:}];
        case 'Hind(A|B)'   , required_distributions = [required_distributions, entropy_distributions.Hind_A_B{:}];
        case 'Chi(A)'      , required_distributions = [required_distributions, entropy_distributions.Chi_A{:}];
        case 'Hsh(A)'      , required_distributions = [required_distributions, entropy_distributions.Hsh_A{:}];
        case 'Hsh(A|B)'    , required_distributions = [required_distributions, entropy_distributions.Hsh_A_B{:}];
        case 'Hnsb(A)'     , required_distributions = [required_distributions, entropy_distributions.Hnsb_A{:}];
        case 'Hnsb(B)'     , required_distributions = [required_distributions, entropy_distributions.Hnsb_B{:}];
        case 'Hnsb(A,B)'   , required_distributions = [required_distributions, entropy_distributions.Hnsb_AB{:}];
    end
end
required_distributions = unique(required_distributions);
prob_dists = prob_estimator(inputs_1d, required_distributions, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                      Step 4.B: Compute requested Entropies                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropies         = cell(1, length(reqOutputs));
entropies_plugin  = cell(1, length(reqOutputs));

symtoolbox = false; % any(strcmp('Symbolic Math Toolbox', {ver().Name}));

for t = 1:nTimepoints
    for i = 1:length(indices)
        idx = indices(i);
        switch possibleOutputs{idx}

            case 'H(A)'
                P_A = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                if symtoolbox
                    P_A = vpa(P_A); P_lin_log = P_A .* log(P_A) / log(vpa(2));
                else
                    P_lin_log = P_A .* log(P_A) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));

                if strcmp(opts.bias, 'bub')
                    bias = bub(nTrials * P_A);
                elseif strcmp(opts.bias, 'pt')
                    bias = (nTimepoints==1) * pt(inputs_1d{1}, nbinsA, nTrials) + ...
                           (nTimepoints>1)  * pt(inputs_1d{1}(1,t,:), nbinsA, nTrials);
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'H(B)'
                P_B = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                if symtoolbox
                    P_B = vpa(P_B); P_lin_log = P_B .* log(P_B) / log(vpa(2));
                else
                    P_lin_log = P_B .* log(P_B) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));

                if strcmp(opts.bias, 'bub')
                    bias = bub(nTrials * P_B);
                elseif strcmp(opts.bias, 'pt')
                    bias = (nTimepoints==1) * pt(inputs_1d{2}, nbinsB, nTrials) + ...
                           (nTimepoints>1)  * pt(inputs_1d{2}(1,t,:), nbinsB, nTrials);
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'H(A|B)'
                P_A_given_B = prob_dists{t, strcmp(required_distributions, 'P(A|B)')};
                P_B         = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                if symtoolbox
                    P_A_given_B = vpa(P_A_given_B); P_B = vpa(P_B);
                    P_lin_log = (P_B.' .* P_A_given_B) .* log(P_A_given_B) / log(vpa(2));
                else
                    P_lin_log = (P_B.' .* P_A_given_B) .* log(P_A_given_B) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));

                if strcmp(opts.bias, 'bub')
                    P_A  = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                    pAB  = P_A_given_B .* P_B.';   % joint
                    bias = bub(nTrials*pAB(:)) - bub(nTrials*P_A(:));
                elseif strcmp(opts.bias, 'pt')
                    bias = 0;
                    if nTimepoints == 1
                        uB = unique(inputs_1d{2}); nb = length(uB);
                        for b_i = 1:nb
                            A_tmp = inputs_1d{1}(inputs_1d{2} == uB(b_i));
                            bias = bias + pt(A_tmp, numel(unique(A_tmp)), nTrials);
                        end
                    else
                        uB = unique(inputs_1d{2}(1,t,:)); nb = length(uB);
                        for b_i = 1:nb
                            mask = (inputs_1d{2}(1,t,:) == uB(b_i));
                            A_tmp = inputs_1d{1}(1,t,mask);
                            bias = bias + pt(A_tmp, numel(unique(A_tmp)), nTrials);
                        end
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'Hlin(A)'
                % Plin(A): per-dimension marginals (rows=dims, cols=bins per dim)
                P_lin = prob_dists{t, strcmp(required_distributions, 'Plin(A)')};
                P_lin_log = P_lin .* log(P_lin) / log(2);
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log,'all'));

                if strcmp(opts.bias, 'bub')
                    bias = 0;
                    for r = 1:size(P_lin,1)
                        pr = P_lin(r, :);
                        bias = bias + bub(nTrials * pr(:));
                    end
                elseif strcmp(opts.bias, 'pt')
                    bias = 0;
                    if nTimepoints == 1
                        for r = 1:size(inputs{1},1)
                            Arow  = inputs{1}(r,:);            % 1 x nTrials
                            nArow = numel(unique(Arow(:)));
                            bias  = bias + pt(Arow, nArow, nTrials);
                        end
                    else
                        for r = 1:size(inputs{1},1)
                            Arow  = squeeze(inputs{1}(r,t,:)); % nTrials x 1
                            nArow = numel(unique(Arow(:)));
                            bias  = bias + pt(Arow, nArow, nTrials);
                        end
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'Hind(A)'
                P_indA = prob_dists{t, strcmp(required_distributions, 'Pind(A)')};
                if symtoolbox
                    P_indA = vpa(P_indA); P_lin_log = P_indA .* log(P_indA) / log(vpa(2));
                else
                    P_lin_log = P_indA .* log(P_indA) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));

                if strcmp(opts.bias, 'pt')
                    bias = 0;
                    if nTimepoints == 1
                        uB = unique(inputs_1d{2}); nb = length(uB);
                        for b_i = 1:nb
                            A_tmp = inputs_1d{1}(inputs_1d{2} == uB(b_i));
                            bias  = bias + pt(A_tmp, nbinsA, nTrials);
                        end
                    else
                        uB = unique(inputs_1d{2}(1,t,:)); nb = length(uB);
                        for b_i = 1:nb
                            mask = (inputs_1d{2}(1,t,:) == uB(b_i));
                            A_tmp = inputs_1d{1}(1,t,mask);
                            bias  = bias + pt(A_tmp, nbinsA, nTrials);
                        end
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'Hind(A|B)'
                P_indAB = prob_dists{t, strcmp(required_distributions, 'Pind(A|B)')};
                P_B     = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                if symtoolbox
                    P_indAB = vpa(P_indAB); P_B = vpa(P_B);
                    P_lin_log = (P_B.' .* P_indAB) .* log(P_indAB) / log(vpa(2));
                else
                    P_lin_log = (P_B.' .* P_indAB) .* log(P_indAB) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));
                bias = 0;
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'Chi(A)'
                P_A    = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                P_indA = prob_dists{t, strcmp(required_distributions, 'Pind(A)')};

                if symtoolbox
                    P_A = vpa(P_A); P_indA = vpa(P_indA);
                    P_lin_log = P_A .* log(P_indA) / log(vpa(2));
                else
                    P_lin_log = P_A .* log(P_indA) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));
                entropies{i}(1,t) = entropies_plugin{i}(1,t); % no bias here

            case 'Hsh(A)'
                P_shA = prob_dists{t, strcmp(required_distributions, 'Psh(A)')};
                if symtoolbox
                    P_shA = vpa(P_shA); P_lin_log = P_shA .* log(P_shA) / log(vpa(2));
                else
                    P_lin_log = P_shA .* log(P_shA) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));

                if strcmp(opts.bias, 'bub')
                    bias = bub(nTrials * P_shA);
                elseif strcmp(opts.bias, 'pt')
                    bias = (nTimepoints==1) * pt(inputs_1d{1}, nbinsA, nTrials) + ...
                           (nTimepoints>1)  * pt(inputs_1d{1}(1,t,:), nbinsA, nTrials);
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'Hsh(A|B)'
                P_shA_given_B = prob_dists{t, strcmp(required_distributions, 'Psh(A|B)')};
                P_B           = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                if symtoolbox
                    P_shA_given_B = vpa(P_shA_given_B); P_B = vpa(P_B);
                    P_lin_log = (P_B.' .* P_shA_given_B) .* log(P_shA_given_B) / log(vpa(2));
                else
                    P_lin_log = (P_B.' .* P_shA_given_B) .* log(P_shA_given_B) / log(2);
                end
                P_lin_log(~isfinite(P_lin_log)) = 0;
                entropies_plugin{i}(1,t) = double(-sum(P_lin_log(:)));

                if strcmp(opts.bias, 'bub')
                    P_shA  = prob_dists{t, strcmp(required_distributions, 'Psh(A)')};
                    PshAB  = P_shA_given_B .* P_B.';  % joint
                    bias   = bub(nTrials*PshAB(:)) - bub(nTrials*P_shA(:));
                elseif strcmp(opts.bias, 'pt')
                    bias = 0;
                    if nTimepoints == 1
                        uB = unique(inputs_1d{2}); nb = length(uB);
                        for b_i = 1:nb
                            A_tmp = inputs_1d{1}(inputs_1d{2} == uB(b_i));
                            bias  = bias + pt(A_tmp, nbinsA, nTrials);
                        end
                    else
                        uB = unique(inputs_1d{2}(1,t,:)); nb = length(uB);
                        for b_i = 1:nb
                            mask = (inputs_1d{2}(1,t,:) == uB(b_i));
                            A_tmp = inputs_1d{1}(1,t,mask);
                            bias  = bias + pt(A_tmp, nbinsA, nTrials);
                        end
                    end
                else
                    bias = 0;
                end
                entropies{i}(1,t) = entropies_plugin{i}(1,t) - bias;

            case 'Hnsb(A)'
                qfun = 1; precision = .1;
                P_A  = prob_dists{t, strcmp(required_distributions, 'P(A)')};
                nb   = P_A' * nTrials; K = length(nb);
                nxa  = nb(nb>0); kxb = ones(size(nxa));
                [Sb_nsb, ~, ~, ~, ~, S_ml, ~] = find_nsb_entropy(kxb, nxa, K, precision, qfun);
                entropies_plugin{i}(1,t) = S_ml;
                entropies{i}(1,t)        = Sb_nsb;

            case 'Hnsb(B)'
                qfun = 1; precision = .1;
                P_B  = prob_dists{t, strcmp(required_distributions, 'P(B)')};
                nb   = P_B' * nTrials; K = length(nb);
                nxb  = nb(nb>0); kxb = ones(size(nxb));
                [Sb_nsb, ~, ~, ~, ~, S_ml, ~] = find_nsb_entropy(kxb, nxb, K, precision, qfun);
                entropies_plugin{i}(1,t) = S_ml;
                entropies{i}(1,t)        = Sb_nsb;

            case 'Hnsb(A,B)'
                qfun = 1; precision = .1;
                P_AB = prob_dists{t, strcmp(required_distributions, 'P(A,B)')};
                nab  = P_AB(:)' * nTrials; K = length(nab);
                nxab = nab(nab>0); kxab = ones(size(nxab));
                [Sab_nsb, ~, ~, ~, ~, S_ml, ~] = find_nsb_entropy(kxab, nxab, K, precision, qfun);
                entropies_plugin{i}(1,t) = S_ml;
                entropies{i}(1,t)        = Sab_nsb;

        end
    end
end
end
