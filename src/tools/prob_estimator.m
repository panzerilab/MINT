function prob_dists = prob_estimator(inputs, reqOutputs, opts)
% *function prob_dists = prob_estimator(inputs, reqOutputs, opts)*
%
% The prob_estimator function calculates various probability distributions based on
% inputs A and B. It computes joint and conditional probabilities, as well
% as independent distributions.
%
% Inputs:
%   - inputs: A cell array containing X data sets:
%             - inputs{1}: First data input (A) with dimensions
%                          nDims [X nTimepoints] X nTrials
%             - inputs{2}: Second data input (B) with dimensions
%                          nDims [X nTimepoints] X nTrials)
%             - inputs{3}: Thirs data input (C) with dimensions
%                          nDims [X nTimepoints] X nTrials)
%                                   ...
%             - inputs{N}: data input N with dimensions
%                          nDims [X nTimepoints] X nTrials)
%
%   - reqOutputs: A cell array of strings specifying which probability distributions
%              to compute. Possible reqOutputs include:
%               - 'P(A)'         : Joint probability distribution of A.
%               - 'P(B)'         : Joint probability distribution of B.
%               - 'P(A,B)'       : Joint probability distribution of A and B.
%               - 'P(A|B)'       : Conditional probability distribution of A given B.
%               - 'Pind(A)'      : Independent joint probability distribution of A.
%               - 'Pind(A|B)'    : Independent conditional probability distribution of A given B.
%               - 'Psh(A|B)'     : Shuffled conditional probability distribution of A given B.
%               - 'Psh(A)'       : Shuffled probability distribution of A.
%               - 'Plin(A)'      : Independent joint probability distribution of A computed linearly.
%
%               - 'P(A,B,C)'     : Multidim Joint probability distribution of A, B and C.
%               - 'P(all)'       : Multidim Joint probability distribution of all input vars.   
%
% Outputs:
%   - prob_dists: A cell array containing the estimated probability distributions
%                 as specified in the reqOutputs argument.
%
% Note:
% Input A and B can represent multiple trials and dimensions concatenated along the
% first dimension. This allows the analysis of interactions between different signals
% and their influence on the probability distributions being studied.
%
% EXAMPLE
% Suppose we have two time series data sets A and B with the following structures:
% A = randn( nDims, nTimepoints, nTrials);% Random data for time series A
% B = randn( nDims, nTimepoints, nTrials);% Random data for time series B
%
% To compute the probability distributions, the function can be called as:
% prob_dists = prob_estimator({A, B}, {'P(A)', 'P(B)', 'P(A|B)'}, opts);
%
% Here, 'opts' represents additional options you may want to include, such as
% specifying the delay (tau), present timepoint (tpres), number of bins (n_bins),
% and other parameters as needed.

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

%% --- Which outputs need linear/independent pieces?
if (length(inputs)>2 && ~any(strcmp(reqOutputs,'P(A,B,C)'))) || (length(inputs)>3 && any(strcmp(reqOutputs,'P(A,B,C)')))
    A_nd = inputs{end};
else
    A_nd = inputs{1};
end
if any(strcmp(reqOutputs,'P(A,B,C)')) || any(strcmp(reqOutputs,'Plin(A)')) || any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Pind(A)')) || any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Psh(A)'))
    % if (length(inputs)>2 && ~any(strcmp(reqOutputs,'P(A,B,C)'))) || (length(inputs)>3 && any(strcmp(reqOutputs,'P(A,B,C)')))
    %     A_nd = inputs{end};
    % else
    %     A_nd = inputs{1};
    % end
    needLin = true;
else
    needLin = false;
end

%% --- Load & normalize inputs (1-based codes), reduce dims to 1D indices for convenience
prob_dists = cell(size(reqOutputs));
nTimepoints = 0;

for var = 1:length(inputs)
    letter = char(64 + var);
    data.(letter) = inputs{var};

    % 1-based integer coding
    min_val = min(data.(letter), [], 'all');
    if min_val ~= 1
        data.(letter) = data.(letter) - min_val + 1;
    end

    if ndims(data.(letter)) == 2
        [nDim, nTrials] = size(data.(letter));
        nTimepointsCur = 1;
    elseif ndims(data.(letter)) == 3
        [nDim, nTimepointsCur, nTrials] = size(data.(letter));
    else
        error(['probEstim: Input ' letter ' is in the wrong shape.']);
    end

    nTimepoints = max(nTimepoints, nTimepointsCur);
    nTimepoints_data.(letter) = nTimepointsCur;

    % If multidim along rows, keep both: the full multi-dim values & the 1D index view.
    % if nDim > 1
    %     data_full.(letter) = data.(letter);  % trials x time x dims
    %     [data_1d.(letter), resps_grid]   = reduce_dim_for_time(data.(letter));% 1 x time x trials (indices)
    % else
    %     data_full.(letter) = permute(data.(letter), [2 3 1]);  % trials x time x 1
    %     data_1d.(letter)   = data.(letter);                    % already 1D
    % end
    if nDim > 1
        % dims x T x N  ->  N x T x dims
        data_full.(letter) = permute(data.(letter), [3 2 1]);
        data_1d.(letter)   = reduce_dim_for_time(data.(letter)); % update this helper (see #3 below)
    else
        % 1 x N  or 1 x T x N  ->  N x T x 1
        data_full.(letter) = permute(data.(letter), [3 2 1]);
        data_1d.(letter)   = squeeze(permute(data.(letter), [3 2 1])); % N x T
    end

end

%% --- Main loop over time
for t = 1:max(1, nTimepoints)
    clear A_t B_t C_t FullA_t p_A p_B p_AB p_A_B psh_A psh_A_B p_all p_ABC ...
          PindA PindAB PinabA PinabAB PlinA;

    % Extract per-variable slices at this time
    for var = 1:length(inputs)
        letter = char(64 + var);

        if nTimepoints_data.(letter) > 1
            % trials x dims (full), trials x 1 (1d)
            Full_t  = data_full.(letter)(:, t, :);     % nTrials x nDimsVar
            Index_t = data_1d.(letter)(:, t);          % nTrials x 1
        else
            Full_t  = squeeze(data_full.(letter));
            if size(Full_t,1) ~= nTrials, Full_t = Full_t.'; end
            if ndims(inputs{var}) > 2
                Index_t = squeeze(data_1d.(letter)(:, 1));
            else
                Index_t = data_1d.(letter).';
                if size(Index_t,1) ~= nTrials, Index_t = Index_t.'; end
            end
        end

        if size(Full_t,1) ~= nTrials, Full_t = Full_t.'; end
        if size(Index_t,1) ~= nTrials, Index_t = Index_t.'; end

        data_full_t.(letter) = Full_t;   % nTrials x nDimsVar
        data_1d_t.(letter)   = Index_t;  % nTrials x 1
    end

    % For models that need the full multi-dim A
    if needLin
        % Decide whether "A" should be taken from the last input or from A
        useLastAsA = (length(inputs) > 2 && ~any(strcmp(reqOutputs,'P(A,B,C)'))) || ...
                     (length(inputs) > 3 &&  any(strcmp(reqOutputs,'P(A,B,C)')));
    
        if useLastAsA
            lastLetter = char(64 + length(inputs));   % 'A','B','C',... for the last input
            FullA_t = A_nd; %data_full_t.(lastLetter);       % nTrials x nDims_of_last_input
        else
            FullA_t = A_nd; %data_full_t.A;                  % nTrials x nDims_of_A
        end
    end

    % Convenience local A/B/C matrices (nTrials x dims)
    % A_full = data_full_t.C;                           % nTrials x KA
    if (isfield(data_full_t,'C') & needLin), A_full = data_full_t.C; else, A_full = data_full_t.A; end
    if isfield(data_full_t,'B'), B_full = data_full_t.B; else, B_full = []; end
    if isfield(data_full_t,'C'), C_full = data_full_t.C; else, C_full = []; end

    % ----- Build full supports (Cartesian) we will reuse -----
    % A support (all patterns over A's dims)

    % B support (if present) uses 1D values only (if B has multiple dims, we also go full Cartesian)
    if ~isempty(B_full)
        [B_patterns, B_levels] = build_support(B_full); % M_B x KB
        % In the usual 1D-B case, M_B == number of unique(B)
    end
    if ~isempty(C_full)
        [C_patterns, C_levels] = build_support(C_full);
    end
    
    if needLin
        [A_patterns, A_levels] = build_support(A_full);   % M_A x KA ; A_levels is 1xKA cell of unique values
        A_idx = map_rows_to_patterns(A_full, A_patterns);     % nTrials x 1 in 1..M_A
    else
        A_patterns = (1:max(data_1d_t.A))';
        A_levels = data_1d_t.A;
        A_idx = A_levels;
    end
    % Precompute row -> pattern indices for speed
    if ~isempty(B_full), B_idx = map_rows_to_patterns(B_full, B_patterns); end
    if ~isempty(C_full), C_idx = map_rows_to_patterns(C_full, C_patterns); end

    nTrials_local = size(A_full,1);
    % end

    %% ---------- P(all): joint over ALL inputs (full support) ----------
    if any(strcmp(reqOutputs, 'P(all)'))
        % Build support levels freshly for each input (never rely on A_levels,
        % which is a plain numeric vector when needLin=false).
        [~, pall_A_levels] = build_support(A_full);
        ALL_full = A_full;
        sup_all  = pall_A_levels;
        if ~isempty(B_full)
            [~, pall_B_levels] = build_support(B_full);
            ALL_full = [ALL_full, B_full];           %#ok<AGROW>
            sup_all  = [sup_all, pall_B_levels];     %#ok<AGROW>
        end
        for var = 3:length(inputs)
            letter = char(64 + var);
            [~, Xlev] = build_support(data_full_t.(letter));
            ALL_full = [ALL_full, data_full_t.(letter)];  %#ok<AGROW>
            sup_all  = [sup_all, Xlev];                   %#ok<AGROW>
        end
        [ALL_patterns, ~] = build_support_from_levels(sup_all);
        ALL_idx = map_rows_to_patterns(ALL_full, ALL_patterns);

        p_all = accumarray(ALL_idx, 1, [size(ALL_patterns,1), 1]);
        p_all = p_all / sum(p_all);
        nLevels = cellfun(@numel, sup_all);
        p_all = reshape(p_all, nLevels);
    end

    %% ---------- P(A,B,C) on full support ----------
    if any(strcmp(reqOutputs,'P(A,B,C)'))
        if isempty(B_full) || isempty(C_full)
            error('P(A,B,C) requested but B or C is missing.');
        end
        [A,~,ia] = unique(data_1d_t.A);
        [B,~,ib] = unique(data_1d_t.B);
        [C,~,ic] = unique(data_1d_t.C);
        
        P = accumarray([ia ib ic], 1, [length(A) length(B) length(C)]);
        p_ABC = P / sum(P(:));
    end

    %% ---------- P(B) (full support over B) ----------
    if any(strcmp(reqOutputs,'P(B)')) || any(strcmp(reqOutputs,'P(A|B)')) || any(strcmp(reqOutputs,'P(A,B)')) ...
                                      || any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Pind(A|B)')) ...
                                      || any(strcmp(reqOutputs,'Psh(A)'))
        if isempty(B_full), error('B requested but missing.'); end
        p_B = accumarray(B_idx, 1, [size(B_patterns,1), 1]) / nTrials_local;
    end

    %% Independent per-dimension probability distributions of A (plin_A)
    % Each row k is P(A_k = i), i over that row's own support.
    if any(strcmp(reqOutputs,'Plin(A)'))
        KA = size(A_nd, 1);                        % #dims in A at this time
        nbinsA = zeros(1, KA);
        for k = 1:KA
            nbinsA(k) = numel(unique(A_nd(k,:)));  % support size for dim k
        end
        maxbins = max(nbinsA);
        plin_A  = zeros(KA, maxbins);
    
        for k = 1:KA
            lvls = unique(A_nd(k,:));              % sorted support for dim k
            vals = A_nd(k,:);                      % 1 x nTrials
    
            % Map values to indices on this dim's support
            [~, pos] = ismember(vals, lvls);       % pos in 1..numel(lvls) or 0 if not found
            pos = pos(:);                          % <-- ensure column
            pos = pos(pos > 0);                    % <-- drop any zeros (robustness)
    
            % Counts aligned to lvls
            cnt = accumarray(pos, 1, [numel(lvls), 1]);
    
            % Normalize to probabilities
            pk  = cnt / numel(vals);               % or / nTrials_local if you prefer
            plin_A(k, 1:numel(lvls)) = pk.';       % row k = P(A_k = i)
        end
    end


  
    %% ---------- P(A) (full support over A) ----------
    if any(strcmp(reqOutputs,'P(A)')) || any(strcmp(reqOutputs,'Pind(A)'))
        p_A = accumarray(A_idx, 1, [size(A_patterns,1), 1]) / nTrials_local;
    end

    %% ---------- P(A,B) (full support: |A| x |B|) ----------
    if any(strcmp(reqOutputs,'P(A,B)')) || any(strcmp(reqOutputs,'P(A|B)')) || ...
       any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Psh(A)'))
        if isempty(B_full), error('P(A,B) requested but B is missing.'); end
        lin = sub2ind([size(A_patterns,1), size(B_patterns,1)], A_idx, B_idx);
        countsAB = accumarray(lin, 1, [size(A_patterns,1)*size(B_patterns,1), 1]);
        p_AB = reshape(countsAB, size(A_patterns,1), size(B_patterns,1));
        p_AB = p_AB / sum(p_AB(:));
    end

    %% ---------- P(A|B) (column-normalized P(A,B)) ----------
    if any(strcmp(reqOutputs,'P(A|B)'))
        p_A_B = normalize_cols(p_AB);
    end

    %% ---------- Shuffled models (keep full A support) ----------
    if any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Psh(A)')) %|| any(strcmp(reqOutputs,'Pind(A|B)'))
        % Shuffle A *within* B (preserving B's distribution)
        shuffled_A = shuffle_core(data_1d_t.B, A_full, 1, [1 0]);  % returns nTrials x KA values
        % Map shuffled rows to A's full support
        A_idx_sh = map_rows_to_patterns(shuffled_A, A_patterns);
        % Joint with B on full grid
        lin_sh = sub2ind([size(A_patterns,1), size(B_patterns,1)], A_idx_sh, B_idx);
        counts_sh = accumarray(lin_sh, 1, [size(A_patterns,1)*size(B_patterns,1), 1]);
        psh_AB = reshape(counts_sh, size(A_patterns,1), size(B_patterns,1));
        psh_AB = psh_AB / sum(psh_AB(:));

        psh_A   = sum(psh_AB, 2);
        psh_A_B = normalize_cols(psh_AB);
    end

%% ---------- Independent Pind via implicit-expansion tensors (no kron) ----------
if any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Pind(A)'))
    % Dimensions of A based on current time-slice
    K   = size(A_full, 2);                 % #dims of A
    nAk = cellfun(@numel, A_levels);       % bins per A-dim (aligned to A_levels)
    MA  = prod(nAk);                       % |A| (full Cartesian support)
    Ntr = size(A_full,1);

    % helper to multiply per-dim marginals into a K-D tensor then vectorize


    if ~isempty(B_full)
        % ----- Conditional independent: Pind(A|B=b) = ⊗_k P(A_k|B=b) (without kron) -----
        nB = size(B_patterns, 1);
        PindA_B = zeros(MA, nB);

        for bi = 1:nB
            sel = (B_idx == bi);
            tot = sum(sel);
            if tot == 0, continue; end

            pk = cell(1, K);
            for k = 1:K
                [~, pos] = ismember(A_full(sel, k), A_levels{k});
                pos = pos(:); pos = pos(pos > 0);
                cnt = accumarray(pos, 1, [nAk(k), 1]);
                pk{k} = cnt / tot;                     % column
            end

            PindA_B(:, bi) = tensor_from_pk(pk, nAk);
        end

        % ----- Unconditional independent as B-marginal of the above -----
        if any(strcmp(reqOutputs,'Pind(A)'))
            if ~exist('p_B','var')
                p_B = accumarray(B_idx, 1, [size(B_patterns,1), 1]) / Ntr;
            end
            PindA = PindA_B * p_B;   % |A| x 1
        end

    else
        % ----- No B provided: Pind(A) = ⊗_k P(A_k) (without kron) -----
        if any(strcmp(reqOutputs,'Pind(A)'))
            pk = cell(1, K);
            for k = 1:K
                [~, pos] = ismember(A_full(:, k), A_levels{k});
                pos = pos(:); pos = pos(pos > 0);
                cnt = accumarray(pos, 1, [nAk(k), 1]);
                pk{k} = cnt / Ntr;                      % column
            end
            PindA = tensor_from_pk(pk, nAk);           % |A| x 1
        end

        if any(strcmp(reqOutputs,'Pind(A|B)'))
            error('Pind(A|B) requested but B is missing.');
        end
    end
end



    %% ---------- Collect outputs for this time ----------
    for outidx = 1:length(reqOutputs)
        name = reqOutputs{outidx};
        switch name
            case 'P(A)',         prob_dist_result = p_A;
            case 'Plin(A)',      prob_dist_result = plin_A;
            case 'P(B)',         prob_dist_result = p_B;
            case 'P(A,B)',       prob_dist_result = p_AB;
            case 'P(A|B)',       prob_dist_result = p_A_B;
            case 'Pind(A)',      prob_dist_result = PindA;
            case 'Pind(A|B)',    prob_dist_result = PindA_B;
            case 'Psh(A|B)',     prob_dist_result = psh_A_B;
            case 'Psh(A)',       prob_dist_result = psh_A;
            case 'P(A,B,C)',     prob_dist_result = p_ABC;
            case 'P(all)',       prob_dist_result = p_all;
            otherwise, error('Unknown reqOutput: %s', name);
        end

        if nTimepoints > 1
            prob_dists{t, outidx} = prob_dist_result;
        else
            prob_dists{outidx} = prob_dist_result;
        end
    end
end

%% ==================== helpers ====================

function [idx1d, resps_grid] = reduce_dim_for_time(R)
    % Wanted shape:
    %   - 2D (nDim x nTrials):  nTrials x 1
    %   - 3D (nDim x T x N):    nTrials x T
    if size(R,1) > 1
        if ndims(R) == 2
            [tmp, resps_grid]   = reduce_dim(R, 1);   % 1 x nTrials
            idx1d = tmp.';              % nTrials x 1
        else
            T     = size(R,2);
            N     = size(R,3);
            idx1d = zeros(N, T);        % nTrials x T
            for tt = 0:T-1
                [tmp, resps_grid] = reduce_dim(R(:,tt+1,:), 1);  % 1 x 1 x N
                idx1d(:, tt+1) = squeeze(tmp).';   % N x 1 -> column
            end
        end
    else
        % Already 1D along rows; make sure trials are rows
        if ndims(R) == 2            % 1 x nTrials
            idx1d = R.';            % nTrials x 1
        else                        % 1 x T x nTrials
            idx1d = squeeze(permute(R, [3 2 1])); % nTrials x T
        end
    end
end


function [patterns, levels] = build_support(X)
    % X: nTrials x K (values are integers starting at 1, but not necessarily contiguous)
    K = size(X,2);
    levels = cell(1,K);
    for k = 1:K
        levels{k} = unique(X(:,k));
    end
    patterns = build_support_from_levels(levels);
end

function [patterns, levels] = build_support_from_levels(levels)
    % levels: 1xK cell, each is a (sorted) vector of unique values
    K = numel(levels);
    [grids{1:K}] = ndgrid(levels{:});
    patterns = reshape(cat(K+1, grids{:}), [], K);  % (#patterns x K)
end

function idx = map_rows_to_patterns(X, patterns)
    % Map each row of X (n x K) to row index in 'patterns' (M x K).
    [tf, idx] = ismember(X, patterns, 'rows'); 
    % If any row wasn't in the Cartesian support (shouldn't happen), set 0
    idx(~tf) = 0;
    if any(~tf)
        warning('Some rows not found in full support; assigning zero index.');
        disp(X);
        disp('patterns')
        disp(patterns)
    end
    % Replace zeros with 1 and we'll zero-weight later if needed
    idx(idx==0) = 1;
end

function M = normalize_cols(M)
    s = sum(M,1);
    s(s==0) = 1;
    M = M ./ s;
end

end

function v = tensor_from_pk(pkCell, nAkLocal)
    % pkCell: 1xK cell, each a marginal (vector)
    % nAkLocal: 1xK bins per dim
    K = numel(pkCell);
    sizes = nAkLocal(:).';
    P = 1;

    for kk = 1:K
        vec = pkCell{kk}(:); % columnize
        if numel(vec) ~= sizes(kk)
            error('tensor_from_pk: length mismatch at dim %d (got %d, expected %d).', ...
                  kk, numel(vec), sizes(kk));
        end
        shp = ones(1, max(K,2)); % ensure at least 2D
        shp(kk) = sizes(kk);
        P = P .* reshape(vec, shp);  % implicit expansion
    end

    v = P(:);
    s = sum(v);
    if s > 0, v = v / s; end
end


