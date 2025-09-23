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

if any(strcmp(reqOutputs,'Plin(A)')) || any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Pind(A)')) || any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Psh(A)'))
    if (length(inputs)>2 && ~any(strcmp(reqOutputs,'P(A,B,C)'))) || (length(inputs)>3 && any(strcmp(reqOutputs,'P(A,B,C)')))
        A_nd = inputs{end};
    else
        A_nd =inputs{1};
    end
    needLin = true;
else
    needLin = false;
end

% Preallocate cell array for the results
prob_dists = cell(size(reqOutputs));
nTimepoints = 0;
for var = 1:length(inputs)
    letter = char(64 + var);
    data.(letter) = inputs{var};
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
    if  nDim > 1
        data_1d.(letter) = reduce_dim(data.(letter), 1);
    else
        data_1d.(letter) = data.(letter);
    end
end
 
% Initialize output for each timepoint if necessary
if nTimepoints > 1
    prob_dists_time = cell(nTimepoints, length(reqOutputs));  % Preallocate for each timepoint
end

% Loop over each timepoint if nTimepoints > 1
for t = 1:max(1, nTimepoints)
    % If nTimepoints > 1, extract the data for the current timepoint
    for var = 1:length(inputs)
        letter = char(64 + var);
        currentData_1d = data_1d.(letter);
        nTimepointsCur =  nTimepoints_data.(letter);
        if nTimepointsCur > 1 || (nTimepointsCur == 1 && ndims(inputs{var}) > 2)
            currentData_t = squeeze(currentData_1d(:, t, :));          
            if needLin
                FullA_t = squeeze(int8(A_nd(:, t, :)));
            end
        else
            currentData_t = currentData_1d;  
            if needLin
                FullA_t = squeeze(int8(A_nd));
            end
        end
        if size(currentData_t, 1) ~= nTrials
            currentData_t = currentData_t';  % Transponiere, wenn notwendig
        end
        data_t.(letter) = currentData_t;       
    end

    if exist('FullA_t', 'var')
        if size(FullA_t,1) ~= nTrials
            FullA_t = FullA_t';
        end
        resps_grid = cell(1, size(FullA_t,2));
    end

    A_t = data_t.A;

    if any(strcmp(reqOutputs, 'P(all)'))
        bins = [];
        all_t = A_t;
        for var = 2:length(inputs)
            letter = char(64 + var);
            currentData_t = data_t.(letter);
            all_t = cat(2,all_t, currentData_t);
        end
        for b = 1:size(all_t, 2)
            k = max(all_t(:, b));
            bins = [bins, k];
        end
        prob_dist_tmp = zeros(bins);
        numtrials = size(all_t, 1);
        for nT = 1:numtrials
            index = num2cell(all_t(nT, :));
            prob_dist_tmp(index{:}) = prob_dist_tmp(index{:}) + 1;
        end
        p_all = prob_dist_tmp / sum(prob_dist_tmp(:));
    end

    if any(strcmp(reqOutputs,'P(A,B,C)'))
        bins = [];
        B_t = data_t.B;
        C_t = data_t.C;
        ABC = cat(2, A_t, B_t, C_t);
        for b = 1:size(ABC, 2)
            k = max(ABC(:,b));%k = length(unique(data(i,:))); %
            bins = [bins, k];
        end
        prob_dist_tmp = zeros(bins);
        numtrials=size(ABC,1);
        for nT = 1:numtrials
            index = num2cell(ABC(nT, :));
            prob_dist_tmp(index{:}) = prob_dist_tmp(index{:})+1;
        end
        p_ABC = prob_dist_tmp/sum(prob_dist_tmp(:));
    end

    
    %% Joint Probability Distribution of B (p_B)
    if any(strcmp(reqOutputs,'P(A,B)')) ||any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Pind(A)')) ||any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Psh(A)'))|| any(strcmp(reqOutputs,'P(B)'))|| any(strcmp(reqOutputs,'P(A|B)'))
        B_t = data_t.B;
        unique_values = unique(B_t);
        ranks = 1:length(unique_values); % Rank values in ascending order
        ranked_B = arrayfun(@(x) ranks(unique_values == x), B_t);
        B_t = ranked_B;

        p_B = prob_estimator_simple(B_t);
    end

    %% Independent Joint Probability Distribution of A (plin_A)
    % Compute independent joint probability as the product of marginal distributions
    if any(strcmp(reqOutputs,'Plin(A)')) || any(strcmp(reqOutputs,'Pind(A)'))
        % Marginal Probability Distributions of A (pmarg_A)
        % plin_A = cell(1, nDimA);  % Cell array to store marginal distributions of each dimension of A
        Ac = mat2cell(FullA_t', ones(1,size(FullA_t',1)), size(FullA_t',2));                 % Split Matrix Into Cells By Row
        [~,edgesall] = histcounts(FullA_t, 'BinMethod','integers');
        [hcell,~] = cellfun(@(X) histcounts(X',edgesall, 'Normalization', 'probability'), Ac, 'Uni',0);   % Do ‘histcounts’ On Each Column
        plin_A = cell2mat(hcell);                                         % Recover Numeric Matrix From Cell Array

        % for dim = 1:nDimA
        %     [~, ~, idx_dim] = unique(A(:, dim));
        %     counts_dim = accumarray(idx_dim, 1);
        %     plin_A(dim,:) = counts_dim / nTrials;
        % end
    end


    %% Joint Probability Distribution of A given B (p_A_B)
    if any(strcmp(reqOutputs,'P(A,B)')) ||any(strcmp(reqOutputs,'P(A|B)')) || any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Psh(A)'))
        % B_t = data_t.B;
        
        p_AB = prob_estimator_simple([A_t, B_t]);
    end

    %% Conditional Joint Probability Distribution of A given B (p_A_B)
    if any(strcmp(reqOutputs,'P(A|B)'))
        %p_A_B = (p_AB'./p_B)';
        p_A_B = p_AB ./ sum(p_AB, 1);
    end

    %% Joint Probability Distribution of A (p_A)
    if any(strcmp(reqOutputs,'P(A)')) || any(strcmp(reqOutputs,'Pind(A)'))
        if exist('p_AB','var') == 1
            p_A = sum(p_AB,2);
        else
            p_A = prob_estimator_simple(A_t);
        end
    end
    %% Independent Conditional Joint Probability Distribution of A (pind_A)
    % if any(strcmp(reqOutputs,'Pind(A)'))
    %     num_shuffles = 10;
    %     pind_A_sum = 0;
    %
    %
    %     for i = 1:num_shuffles
    %         shuffled_A = shuffle_core(B_t, FullA_t, 1, [1 0]);
    %         if size(shuffled_A, 2) > 1
    %             shuffled_A_1d = reduce_dim(shuffled_A', 1);
    %             shuffled_A_1d = shuffled_A_1d';
    %         else
    %             shuffled_A_1d = shuffled_A;
    %         end
    %
    %         pind_A_tmp = prob_estimator_simple(shuffled_A_1d);
    %         if i > 1 && size(pind_A_tmp,1)<size(pind_A_sum,1)
    %             pind_A_tmp = [pind_A_tmp; zeros(size(pind_A_sum,1)-size(pind_A_tmp,1),1)];
    %         end
    %         if i > 1 && size(pind_A_tmp,1)>size(pind_A_sum,1)
    %             pind_A_sum = [pind_A_sum; zeros(size(pind_A_tmp,1)-size(pind_A_sum,1),1)];
    %         end
    %         pind_A_sum = pind_A_sum + pind_A_tmp;
    %     end
    %     pind_A = pind_A_sum ./ num_shuffles;
    % end

    %% Shuffled Joint Probability Distribution of A given B (psh_A_B)
    if any(strcmp(reqOutputs,'Psh(A|B)')) || any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Psh(A)')) ||  any(strcmp(reqOutputs,'Pind(A)'))
        shuffled_A = shuffle_core(B_t, FullA_t, 1, [1 0]);  % Initialize a shuffled version of A
        if  size(shuffled_A,2) > 1
            shuffled_A_1d = reduce_dim(shuffled_A',1);
            shuffled_A_1d = shuffled_A_1d';
        else
            shuffled_A_1d = shuffled_A;
        end
        psh_AB = prob_estimator_simple([shuffled_A_1d, B_t]);
        psh_A_B = (psh_AB'./p_B)';
        psh_A = sum(psh_AB,2);
    end

    % %% Independent Conditional Joint Probability Distribution of A given B (pind_A_B)
    % if any(strcmp(reqOutputs,'Pind(A|B)')) || any(strcmp(reqOutputs,'Pind(A)'))
    %     % num_shuffles = 100;
    %     % psh_AB_sum = 0;  % Initialisierung der Summe für psh_AB
    %     % [a1d, patts] = reduce_dim(FullA_t', 1);
    %     % for k = 1:num_shuffles
    %     %     shuffled_A = shuffle_core(B_t, FullA_t, 1, [1 0]);
    %     %     if size(shuffled_A, 2) > 1
    %     %         [shuffled_A_1d, patts_shuff] = reduce_dim(shuffled_A', 1);
    %     %         if ~all(patts' == patts_shuff')
    %     %             disp('hola')
    %     %         end
    %     %         shuffled_A_1d = shuffled_A_1d';
    %     %     else
    %     %         shuffled_A_1d = shuffled_A;
    %     %     end
    %     %     psh_AB_tmp = prob_estimator_simple([shuffled_A_1d, B_t]);
    %     %     if k > 1 && size(psh_AB_tmp,1)<size(psh_AB_sum,1)
    %     %         psh_AB_tmp = [psh_AB_tmp; zeros(size(psh_AB_sum,1)-size(psh_AB_tmp,1),size(psh_AB_tmp,2))];
    %     %     end
    %     %     if k > 1 && size(psh_AB_tmp,1)>size(psh_AB_sum,1)
    %     %         psh_AB_sum = [psh_AB_sum; zeros(size(psh_AB_tmp,1)-size(psh_AB_sum,1),size(psh_AB_tmp,2))];
    %     %     end
    %     %     psh_AB_sum = psh_AB_sum + psh_AB_tmp;
    %     % end
    %     % pind_ABtest = psh_AB_sum ./ num_shuffles;
    %     % pind_A_B = (pind_ABtest' ./ p_B)';
    %     % pind_A = sum(pind_ABtest, 2);
    % 
    %     % UniqueA = unique(A_t);
    %     % UniqueB = unique(B_t);
    %     % plin_A_B = [];
    %     % [~,edgesall] = histcounts(FullA_t, 'BinMethod','integers');
    %     % [a1d, patts] = reduce_dim(FullA_t', 1);
    %     % nbinsA = {};
    %     % for k=1:size(FullA_t,2)
    %     %     nbinsA{k} =  1:length(unique(FullA_t(:,k)));
    %     % end
    %     % dim_to_collapse = size(FullA_t,2);
    %     % [resps_grid{1:dim_to_collapse}] = ndgrid(nbinsA{:});
    %     % resps_grid = reshape(cat(dim_to_collapse+1, resps_grid{:}), [], dim_to_collapse);
    %     % resps_gridc = mat2cell(resps_grid', ones(1,size(resps_grid',1)), size(resps_grid',2));
    %     % pind_A_B= [];
    %     % for k=1:length(UniqueB)
    %     %     stimA_t = FullA_t(B_t==UniqueB(k),:);
    %     %     Ac = mat2cell(stimA_t', ones(1,size(stimA_t',1)), size(stimA_t',2));                 % Split Matrix Into Cells By Row
    %     %     [hcell,~] = cellfun(@(X) histcounts(X',edgesall, 'Normalization', 'probability'), Ac, 'Uni',0);   % Do ‘histcounts’ On Each Column
    %     %     hcell2 =cellfun(@(X, Y) X(Y), hcell, resps_gridc, 'Uni', 0);
    %     %     condprob = cell2mat(hcell2);
    %     %     pind_Ared_b = prod(condprob,1);
    %     %     pind_A_B = [pind_A_B pind_Ared_b'];
    %     % end
    %     % % pind_A_B = plin_A_B;
    %     % pind_AB = pind_A_B;% * p_B(p_B>0);
    %     % for rowi=1:size(pind_AB,1)
    %     %     pind_AB(rowi,:) = pind_A_B(rowi,:) .* p_B';
    %     % end
    %     % zero_positions = find(p_B == 0);
    %     % if length(p_B) > size(pind_A_B,2)
    %     %     disp('hola')
    %     %     for zp=1:length(zero_positions)
    %     %         % Define the row of zeros
    %     %         new_col = zeros(size(pind_AB, 1),1);  % A col of zeros with the same number of rows as A
    %     %         % Specify the location where you want to insert the new col
    %     %         col_to_insert = zero_positions(zp);
    %     %         % Insert the new col
    %     %         if col_to_insert ==1
    %     %             pind_AB = [new_col, pind_AB];
    %     %         elseif col_to_insert == length(p_B)
    %     %             pind_AB = [pind_AB, new_col];
    %     %         else
    %     %             pind_AB = [pind_AB(:, 1:col_to_insert-1), new_col, pind_AB(:, col_to_insert:end)];
    %     %         end
    %     % 
    %     %     end
    %     % end
    %     % pind_A = sum(pind_AB, 2);
    % 
    % 
    % end

%% --- Independent Conditional (Pind(A|B)) and its B-marginal (Pind(A))
if any(strcmp(reqOutputs,'Pind(A)')) || any(strcmp(reqOutputs,'Pind(A|B)'))
    % Canonical pattern table (values) used by p_A
    [~, A_patterns] = reduce_dim(FullA_t', 1);     % M x K (values, NOT indices)
    M = size(A_patterns,1);
    K = size(A_patterns,2);

    % Global unique value sets per A-dimension (same convention as reduce_dim)
    resps_all = cell(1,K);
    for k = 1:K
        resps_all{k} = unique(FullA_t(:,k));       % sorted values
    end

    % Map pattern VALUES -> INDICES per dimension (so we can index marginals safely)
    A_patterns_idx = zeros(M, K);
    for k = 1:K
        [tf,pos] = ismember(A_patterns(:,k), resps_all{k});
        % if any pattern had a value not present in resps_all (shouldn't happen), mark 0
        pos(~tf) = 0;
        A_patterns_idx(:,k) = pos;
    end

    % Ensure B_t and p_B exist and are 1..nB ranked (as in your earlier code)
    if ~exist('B_t','var') || ~exist('p_B','var')
        B_t_raw = data_t.B;
        uv = unique(B_t_raw);
        ranks = 1:length(uv);
        B_t = arrayfun(@(x) ranks(uv == x), B_t_raw);
        p_B = prob_estimator_simple(B_t);          % column vector length nB
    end
    nB = numel(p_B);

    % -------- Conditional independent Pind(A|B) --------
    pind_A_B = zeros(M, nB);                       % default zeros
    for bi = 1:nB
        sel = (B_t == bi);
        if ~any(sel), continue; end

        % Per-dimension conditionals aligned to resps_all{k}
        marg_b = cell(1,K);
        total_b = sum(sel);
        for k = 1:K
            vals = FullA_t(sel, k);                        % values when B=bi
            [tf,pos] = ismember(vals, resps_all{k});       % -> indices 1..|resps_all{k}|
            cntk = accumarray(pos(tf), 1, [numel(resps_all{k}), 1]);
            marg_b{k} = cntk / total_b;                    % P(A_k=i | B=bi)
        end

        % Product across dims on the SAME pattern rows (using indices)
        col = ones(M,1);
        for r = 1:M
            pr = 1;
            for k = 1:K
                idxk = A_patterns_idx(r,k);                % index into resps_all{k}
                if idxk < 1 || idxk > numel(marg_b{k}) || marg_b{k}(idxk) == 0
                    pr = 0; break;                         % missing/unseen -> zero
                end
                pr = pr * marg_b{k}(idxk);
            end
            col(r) = pr;
        end

        % Normalize column if there is mass
        s = sum(col);
        if s > 0
            col = col / s;
        else
            col = zeros(M,1);
        end
        pind_A_B(:, bi) = col;
    end

    % -------- Define Pind(A) as the B-marginal of Pind(A|B) --------
    pind_A = pind_A_B * p_B;                               % M x 1
end




    %% Collate all results
    for outidx = 1:length(reqOutputs)
        switch reqOutputs{outidx}
            case 'P(A)', prob_dist_result = p_A;
            case 'Plin(A)', prob_dist_result = plin_A;
            case 'P(A,B)', prob_dist_result = p_AB;
            case 'P(A|B)', prob_dist_result = p_A_B;
            case 'P(B)', prob_dist_result = p_B;
            case 'Pind(A)', prob_dist_result = pind_A;
            case 'Pind(A|B)', prob_dist_result = pind_A_B;
            case 'Psh(A|B)', prob_dist_result = psh_A_B;
            case 'Psh(A)', prob_dist_result = psh_A;
            case 'P(A,B,C)', prob_dist_result = p_ABC;
            case 'P(all)', prob_dist_result = p_all;
        end


        % If multiple timepoints, store result per timepoint
        if nTimepoints > 1
            prob_dists{t, outidx} = prob_dist_result;
        else
            prob_dists{outidx} = prob_dist_result;  % If single timepoint, store directly
        end
    end
end


function p = prob_estimator_simple(A)
    warning('off', 'all');

    if size(A, 2) == 2
        % For 2D case
        % Using accumarray to compute joint histogram for two variables
        p = accumarray(A, 1);  % Increment A to avoid zero-indexing issues
        p = p / sum(p(:));  % Normalize to get probabilities
    else
        % For 1D case
        % Using accumarray to compute histogram
        p = accumarray(A, 1);  % Increment A to avoid zero-indexing issues
        p = p / sum(p);  % Normalize to get probabilities
    end
    warning('on', 'all');
end

function products = calculate_products(A)
    % This function calculates all possible products across the rows of the input matrix A.

    % Get the number of rows and columns in the input matrix
    [num_rows, num_cols] = size(A);

    % Initialize cell array to hold each row for ndgrid
    row_cells = cell(1, num_rows);

    % Fill the cell array with each row of the matrix
    for i = 1:num_rows
        row_cells{i} = A(i, :);
    end

    % Use ndgrid to create grids for all combinations
    [grids{1:num_rows}] = ndgrid(row_cells{:});

    % Compute the products across all dimensions
    products = 1;
    for i = 1:num_rows
        products = products .* grids{i};
    end

    % Reshape the output to a proper size
    products = reshape(products, [], num_cols^num_rows);
end
end