function prob_dists = prob_estimator(inputs, reqOutputs, opts)
% *function prob_dists = prob_estimator(inputs, reqOutputs, opts)*
% Fully "full-support" (Cartesian) probability estimation.
%
% All probability arrays are computed over the complete Cartesian product
% of per-dimension unique values (full pattern support). Any unobserved
% pattern gets probability 0 but remains present in the output support.

% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT (GPLv3+).

%% --- Which outputs need linear/independent pieces?
if any(strcmp(reqOutputs,'Plin(A)')) || any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Pind(A)')) || any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Psh(A)'))
    if (length(inputs)>2 && ~any(strcmp(reqOutputs,'P(A,B,C)'))) || (length(inputs)>3 && any(strcmp(reqOutputs,'P(A,B,C)')))
        A_nd = inputs{end};
    else
        A_nd = inputs{1};
    end
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
    if nDim > 1
        data_full.(letter) = data.(letter);  % trials x time x dims
        [data_1d.(letter), resps_grid]   = reduce_dim_for_time(data.(letter));% 1 x time x trials (indices)
    else
        data_full.(letter) = permute(data.(letter), [2 3 1]);  % trials x time x 1
        data_1d.(letter)   = data.(letter);                    % already 1D
    end
end

%% --- Multi-time handling
if nTimepoints > 1
    prob_dists_time = cell(nTimepoints, length(reqOutputs));
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
            Full_t  = squeeze(data_full.(letter)(:, t, :));     % nTrials x nDimsVar
            Index_t = squeeze(data_1d.(letter)(:, t));          % nTrials x 1
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
    A_full = data_full_t.C;                           % nTrials x KA
    if isfield(data_full_t,'B'), B_full = data_full_t.B; else, B_full = []; end
    if isfield(data_full_t,'C'), C_full = data_full_t.C; else, C_full = []; end

    % ----- Build full supports (Cartesian) we will reuse -----
    % A support (all patterns over A's dims)
    [A_patterns, A_levels] = build_support(A_nd');   % M_A x KA ; A_levels is 1xKA cell of unique values
    % B support (if present) uses 1D values only (if B has multiple dims, we also go full Cartesian)
    if ~isempty(B_full)
        [B_patterns, B_levels] = build_support(B_full); % M_B x KB
        % In the usual 1D-B case, M_B == number of unique(B)
    end
    if ~isempty(C_full)
        [C_patterns, C_levels] = build_support(C_full);
    end

    % Precompute row -> pattern indices for speed
    A_idx = map_rows_to_patterns(A_nd', A_patterns);     % nTrials x 1 in 1..M_A
    if ~isempty(B_full), B_idx = map_rows_to_patterns(B_full, B_patterns); end
    if ~isempty(C_full), C_idx = map_rows_to_patterns(C_full, C_patterns); end

    nTrials_local = size(A_full,1);

    %% ---------- P(all): joint over ALL inputs (full support) ----------
    if any(strcmp(reqOutputs, 'P(all)'))
        % Concatenate all variables' full values horizontally
        ALL_full = A_full;
        sup_all  = A_levels;
        if ~isempty(B_full), ALL_full = [ALL_full, B_full]; sup_all = [sup_all, B_levels]; end %#ok<AGROW>
        for var = 3:length(inputs)
            letter = char(64 + var);
            [Xpat, Xlev] = build_support(data_full_t.(letter)); %#ok<ASGLU>
            ALL_full = [ALL_full, data_full_t.(letter)];        %#ok<AGROW>
            sup_all  = [sup_all, Xlev];                          %#ok<AGROW>
        end
        [ALL_patterns, ~] = build_support_from_levels(sup_all);
        ALL_idx = map_rows_to_patterns(ALL_full, ALL_patterns);

        p_all = accumarray(ALL_idx, 1, [size(ALL_patterns,1), 1]);
        p_all = p_all / sum(p_all);
    end

    %% ---------- P(A,B,C) on full support ----------
    if any(strcmp(reqOutputs,'P(A,B,C)'))
        if isempty(B_full) || isempty(C_full)
            error('P(A,B,C) requested but B or C is missing.');
        end
        % Build combined support as Cartesian of (A_patterns × B_patterns × C_patterns)
        [ABC_patterns, ~] = cartesian_join({A_patterns, B_patterns, C_patterns});
        ABC_idx = sub2ind([size(A_patterns,1), size(B_patterns,1), size(C_patterns,1)], ...
                          A_idx, B_idx, C_idx);
        counts = accumarray(ABC_idx, 1, [numel(ABC_patterns(:,1))*0 + prod([size(A_patterns,1), size(B_patterns,1), size(C_patterns,1)]), 1]);
        % Reshape back to 3D then to vector (keep order A x B x C)
        p_ABC = reshape(counts, size(A_patterns,1), size(B_patterns,1), size(C_patterns,1));
        p_ABC = p_ABC / sum(p_ABC(:));
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
    [tf, idx] = ismember(X, patterns, 'rows'); %#ok<ASGLU>
    % If any row wasn't in the Cartesian support (shouldn't happen), set 0
    idx(~tf) = 0;
    if any(~tf)
        warning('Some rows not found in full support; assigning zero index.');
    end
    % Replace zeros with 1 and we'll zero-weight later if needed
    idx(idx==0) = 1;
end

function M = normalize_cols(M)
    s = sum(M,1);
    s(s==0) = 1;
    M = M ./ s;
end

function v = normvec(v)
    s = sum(v);
    if s > 0, v = v / s; end
end

end
        % Helper: balanced Kronecker product over a cell array of column vectors
function v = kronall(vecs)
    if numel(vecs) == 1
        v = vecs{1};
    else
        mid = floor(numel(vecs) / 2);
        v = kron(kronall(vecs(1:mid)), kronall(vecs(mid+1:end)));
    end
end

function v = tensor_from_pk(pkCell, nAkLocal)
    P = 1;                                     % starts scalar, expands
    for kk = 1:numel(pkCell)
        shp       = ones(1, numel(pkCell));
        shp(kk)   = nAkLocal(kk);
        P         = P .* reshape(pkCell{kk}, shp);  % implicit expansion
    end
    v = P(:);
    s = sum(v); if s > 0, v = v / s; end
end
