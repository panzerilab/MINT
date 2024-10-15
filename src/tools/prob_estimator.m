function prob_dists = prob_estimator(inputs, outputs, opts)
%%% *function prob_dists = prob_estimator(inputs, outputs, opts)*
%%%
%%% The prob_estimator function calculates various probability distributions based on 
%%% inputs A and B. It computes joint and conditional probabilities, as well 
%%% as independent distributions.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing two time series data sets:
%%%             - inputs{1}: First data input (A) with dimensions 
%%%                          nDims [X nTimepoints] X nTrials 
%%%             - inputs{2}: Second data input (B) with dimensions 
%%%                          nDims [X nTimepoints] X nTrials)
%%%
%%%   - outputs: A cell array of strings specifying which probability distributions 
%%%              to compute. Possible outputs include:
%%%               - 'P(A)'         : Joint probability distribution of A.
%%%               - 'P(B)'         : Joint probability distribution of B.
%%%               - 'P(A,B)'       : Joint probability distribution of A and B.
%%%               - 'P(A|B)'       : Conditional probability distribution of A given B.
%%%               - 'Pind(A)'      : Independent joint probability distribution of A.
%%%               - 'Pind(A|B)'    : Independent conditional probability distribution of A given B.
%%%               - 'Psh(A|B)'     : Shuffled conditional probability distribution of A given B.
%%%               - 'Psh(A)'        : Shuffled probability distribution of A.     
%%%               - 'Plin(A)'      : Independent joint probability distribution of A computed linearly.
%%%
%%%               - 'P(A,B,C)'     : Multidim Joint probability distribution of A, B and C.
%%%
%%% Outputs:
%%%   - prob_dists: A cell array containing the estimated probability distributions 
%%%                 as specified in the outputs argument.
%%%
%%% Note: 
%%% Input A and B can represent multiple trials and dimensions concatenated along the 
%%% first dimension. This allows the analysis of interactions between different signals 
%%% and their influence on the probability distributions being studied.
%%%
%%% EXAMPLE
%%% Suppose we have two time series data sets A and B with the following structures:
%%% A = randn( nDims, nTimepoints, nTrials);% Random data for time series A
%%% B = randn( nDims, nTimepoints, nTrials);% Random data for time series B
%%%
%%% To compute the probability distributions, the function can be called as:
%%% prob_dists = prob_estimator({A, B}, {'P(A)', 'P(B)', 'P(A|B)'}, opts);
%%%
%%% Here, 'opts' represents additional options you may want to include, such as 
%%% specifying the delay (tau), present timepoint (tpres), number of bins (n_bins), 
%%% and other parameters as needed.

% possibleOutputs = {'P(A)', 'P(B)', 'P(A,B)', 'P(A|B)', 'Pind(A)', 'Pind(A|B)','Psh(A|B)', 'Psh(A)', 'Plin(A)', 'P(A,B,C)'};
% [isMember, indices] = ismember(outputs, possibleOutputs);
% if any(~isMember)
%     nonMembers = outputs(~isMember);
%     msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
%     error('probEstimator:invalidOutput', msg);
% end


if any(strcmp(outputs,'Plin(A)')) 
    if (length(inputs)>2 && ~any(strcmp(outputs,'P(A,B,C)'))) || (length(inputs)>3 && any(strcmp(outputs,'P(A,B,C)')))
        A_nd = inputs{end}';
    else
        A_nd =inputs{1}';
    end
    needLin = true;
else
    needLin = false;
end

% Transpose A and B for easier manipulation (nTrials x nDimA/nDimS)
A = inputs{1};
B = inputs{2};

% Preallocate cell array for the results
prob_dists = cell(size(outputs));

%check for zeros 
min_valA = min(A,[],'all');
min_valB = min(B,[],'all');
if min_valA ~=1
    A = A - min_valA + 1;
end 
if min_valB ~=1
    B = B - min_valB + 1;
end 

% Get the sizes
if ndims(A) == 2
    A = A';
    B = B';
    [nTrials, nDimA] = size(A);
    [~, nDimB] = size(B);
    nTimepoints = 1;
elseif ndims(A) == 3
    [nTrials, nTimepoints, nDimA] = size(A);
    [~, ~, nDimB] = size(B);
else
    error('probEstim: A is in the wrong shape.')
end
if ~isfield(opts, 'multidim')
    opts.multidim = false;
end 

% Initialize output for each timepoint if necessary
if nTimepoints > 1
    prob_dists_time = cell(nTimepoints, length(outputs));  % Preallocate for each timepoint
end

if nTimepoints==1 && any(strcmp(outputs,'P(A,B,C)'))
    bins = [];
    C = inputs{3}';
    min_valC = min(C,[],'all');
    if min_valC ~=1
        C = C - min_valC + 1;
    end
    dimsC = size(C);
    if dimsC(1) ~= nTrials
        error('prob_estim: C is in the wrong shape.')
    end 
    ABC = cat(2, A, B, C);

    for i = 1:size(ABC, 2)
        k = max(ABC(:,i));%k = length(unique(data(i,:))); %
        bins = [bins, k];
    end
    prob_dist_tmp = zeros(bins);
    numtrials=size(ABC,1);
    for i = 1:numtrials
        index = num2cell(ABC(i, :));
        prob_dist_tmp(index{:}) = prob_dist_tmp(index{:})+1;
    end
    prob_dists = prob_dist_tmp/sum(prob_dist_tmp(:));
    return
end

if  nDimB > 1
    B_1d = (reduce_dim(B',1))';
else
    B_1d = B;
end

if  nDimA > 1
    A_1d = (reduce_dim(A',1))';
else
    A_1d = A;
end

% Loop over each timepoint if nTimepoints > 1
for t = 1:max(1, nTimepoints)
    % If nTimepoints > 1, extract the data for the current timepoint
    if nTimepoints > 1
        A_t = squeeze(A_1d(:, t, :));  % Get A for current timepoint (nTrials x nDimA)
        B_t = squeeze(B_1d(:, t, :));  % Get B for current timepoint (nTrials x nDimB)
        if needLin
            FullA_t = squeeze(A_nd(:, t, :));
        end
    else
        A_t = A_1d;  % Use the whole data if nTimepoints == 1
        B_t = B_1d;
        if  needLin
            FullA_t = A_nd;
        end
    end


    %% Joint Probability Distribution of B (p_B)
    if any(strcmp(outputs,'Pind(A|B)')) || any(strcmp(outputs,'Psh(A|B)')) || any(strcmp(outputs,'Psh(A)'))|| any(strcmp(outputs,'P(B)'))|| any(strcmp(outputs,'P(A|B)'))
        p_B = prob_estimator_simple(B_t);
    end

    %% Independent Joint Probability Distribution of A (plin_A)
    % Compute independent joint probability as the product of marginal distributions
    if any(strcmp(outputs,'Plin(A)')) || any(strcmp(outputs,'Pind(A)'))
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
    if any(strcmp(outputs,'P(A,B)')) ||any(strcmp(outputs,'P(A|B)')) || any(strcmp(outputs,'Psh(A|B)')) || any(strcmp(outputs,'Pind(A|B)')) || any(strcmp(outputs,'Psh(A)'))
        p_AB = prob_estimator_simple([A, B_t]);
    end

    %% Conditional Joint Probability Distribution of A given B (p_A_B)
    if any(strcmp(outputs,'P(A|B)'))
        %p_A_B = (p_AB'./p_B)';
        p_A_B = p_AB ./ sum(p_AB, 1);
    end

    %% Joint Probability Distribution of A (p_A)
    if any(strcmp(outputs,'P(A)')) || any(strcmp(outputs,'Pind(A)'))
        if exist('p_AB','var') == 1
            p_A = sum(p_AB,2);
        else
            p_A = prob_estimator_simple(A_1d);
        end
    end
    %% Independent Conditional Joint Probability Distribution of A (pind_A)
    % if any(strcmp(outputs,'Pind(A)'))
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
    if any(strcmp(outputs,'Psh(A|B)')) || any(strcmp(outputs,'Pind(A|B)')) || any(strcmp(outputs,'Psh(A)'))
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

        % 
        % num_shuffles = 1;
        % psh_AB_sum = 0;  % Initialisierung der Summe für psh_AB
        % for i = 1:num_shuffles
        %     shuffled_A = shuffle_core(B_t, FullA_t, 1, [1 0]);  
        %     if size(shuffled_A, 2) > 1
        %         shuffled_A_1d = reduce_dim(shuffled_A', 1);  
        %         shuffled_A_1d = shuffled_A_1d';
        %     else
        %         shuffled_A_1d = shuffled_A;
        %     end
        %     psh_AB_tmp = prob_estimator_simple([shuffled_A_1d, B_t]);
        %     if i > 1 && size(psh_AB_tmp,1)<size(psh_AB_sum,1)
        %         psh_AB_tmp = [psh_AB_tmp; zeros(size(psh_AB_sum,1)-size(psh_AB_tmp,1),size(psh_AB_tmp,2))];
        %     end
        %     if i > 1 && size(psh_AB_tmp,1)>size(psh_AB_sum,1)
        %         psh_AB_sum = [psh_AB_sum; zeros(size(psh_AB_tmp,1)-size(psh_AB_sum,1),size(psh_AB_tmp,2))];
        %     end
        %     psh_AB_sum = psh_AB_sum + psh_AB_tmp;
        % end
        % psh_AB = psh_AB_sum ./ num_shuffles;
        % psh_A_B = (psh_AB' ./ p_B)';  % p(A|B)
        % psh_A = sum(psh_AB, 2);
    end

    %% Independent Conditional Joint Probability Distribution of A given B (pind_A_B)
    if any(strcmp(outputs,'Pind(A|B)')) || any(strcmp(outputs,'Pind(A)'))
        num_shuffles = 30;
        psh_AB_sum = 0;  % Initialisierung der Summe für psh_AB
        for i = 1:num_shuffles
            shuffled_A = shuffle_core(B_t, FullA_t, 1, [1 0]);  
            if size(shuffled_A, 2) > 1
                shuffled_A_1d = reduce_dim(shuffled_A', 1);  
                shuffled_A_1d = shuffled_A_1d';
            else
                shuffled_A_1d = shuffled_A;
            end
            psh_AB_tmp = prob_estimator_simple([shuffled_A_1d, B_t]);
            if i > 1 && size(psh_AB_tmp,1)<size(psh_AB_sum,1)
                psh_AB_tmp = [psh_AB_tmp; zeros(size(psh_AB_sum,1)-size(psh_AB_tmp,1),size(psh_AB_tmp,2))];
            end
            if i > 1 && size(psh_AB_tmp,1)>size(psh_AB_sum,1)
                psh_AB_sum = [psh_AB_sum; zeros(size(psh_AB_tmp,1)-size(psh_AB_sum,1),size(psh_AB_tmp,2))];
            end
            psh_AB_sum = psh_AB_sum + psh_AB_tmp;
        end
        pind_AB = psh_AB_sum ./ num_shuffles;
        pind_A_Btest = (pind_AB' ./ p_B)';
        pind_Atest = sum(pind_AB, 2);

        UniqueA = unique(A_t);
        UniqueB = unique(B_t);
        plin_A_B = [];
        [~,edgesall] = histcounts(FullA_t, 'BinMethod','integers');

        nbinsA = {};
        for i=1:size(FullA_t,2)
            nbinsA{i} =  1:length(unique(FullA_t(:,i)));
        end 
        dim_to_collapse = size(FullA_t,2);
        [resps_grid{1:dim_to_collapse}] = ndgrid(nbinsA{:});
        resps_grid = reshape(cat(dim_to_collapse+1, resps_grid{:}), [], dim_to_collapse);
        resps_gridc = mat2cell(resps_grid', ones(1,size(resps_grid',1)), size(resps_grid',2));
        pind_A_B= [];
        for i=1:length(UniqueB)
            stimA_t = FullA_t(B_t==UniqueB(i),:);
            Ac = mat2cell(stimA_t', ones(1,size(stimA_t',1)), size(stimA_t',2));                 % Split Matrix Into Cells By Row
            [hcell,~] = cellfun(@(X) histcounts(X',edgesall, 'Normalization', 'probability'), Ac, 'Uni',0);   % Do ‘histcounts’ On Each Column
            hcell2 =cellfun(@(X, Y) X(Y), hcell, resps_gridc, 'Uni', 0);
            condprob = cell2mat(hcell2);
            pind_Ared_b = prod(condprob,1);
            pind_A_B = [pind_A_B pind_Ared_b'];
        end
        % pind_A_B = plin_A_B;
        pind_AB = pind_A_B * p_B;
        pind_A = sum(pind_AB, 2);
    end



    %% Collate all results
    for i = 1:length(outputs)
        switch outputs{i}
            case 'P(A)', prob_dist_result = p_A;
            case 'Plin(A)', prob_dist_result = plin_A;
            case 'P(A,B)', prob_dist_result = p_AB;
            case 'P(A|B)', prob_dist_result = p_A_B;
            case 'P(B)', prob_dist_result = p_B;
            case 'Pind(A)', prob_dist_result = pind_A;
            case 'Pind(A|B)', prob_dist_result = pind_A_B;
            case 'Psh(A|B)', prob_dist_result = psh_A_B;
            case 'Psh(A)', prob_dist_result = psh_A;
        end


        % If multiple timepoints, store result per timepoint
        if nTimepoints > 1
            prob_dists_time{t, i} = prob_dist_result;
        else
            prob_dists{i} = prob_dist_result;  % If single timepoint, store directly
        end
    end
end
end
%
% function p = prob_estimator_simple(A)
%     warning('off','all')
%     if size(A,2)==2
%        p = histcounts2(A(:,1)', A(:,2)', 'Normalization','probability', 'BinMethod','integers');
%     else
%         p = histcounts(A, 'Normalization','probability', 'BinMethod','integers');
%     end
%     warning('on','all')
% end

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