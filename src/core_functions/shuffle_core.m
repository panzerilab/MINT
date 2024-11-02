function shuffled_data = shuffle_core(behav_data, neural_data, consistency, index)
% Shuffle function for neural data with optional consistency constraints.
% This function shuffles neural data based on specified
% conditions, allowing customization through input parameters.

% Input:
%   - behav_data: Behavioral data used for consistency constraints
%   - neural_data: nTrials x nObjects x ..., neural data to be shuffled
%   - consistency: Binary flag (0 or 1) indicating whether to apply
%                  consistency constraints for shuffling across trials (0 for no constraints)
%   - index: Binary vector indicating which dimensions to shuffle

% Output:
%   - shuffled_data: Shuffled neural data based on specified conditions

% Description:
% The function shuffles neural data based on the provided input parameters.
% - If consistency is set to 1, the shuffling is consistent across trials,
%   ensuring that unique values in the behavioral data remain grouped together.
% - The index vector indicates which dimensions should be shuffled. 
%   index(2) == 1 signifies shuffling over the second dimension (neuralObjects),
%   potentially altering the probability distributions of individual objects.

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

if nargin < 4
    msg = "not enough input arguments.";
    error('shuffle_core:notEnoughInput', msg);
end 

nObjects = size(neural_data, 2);
nTrials = size(neural_data, 1);
nDims = ndims(neural_data);
shuffled_data = neural_data;

if index(1)==1
    % check if shuffling across trials should be consistent 
    if consistency == 0
        for col = 1:nObjects
            temp_data = neural_data(:, col, :);
            temp_data = temp_data(randperm(nTrials), :);
            shuffled_data(:, col, :) = temp_data;
        end
    else
        % Consistency constraints based on behavioral data
        [unique_values, ~, ~] = unique(behav_data, 'rows', 'stable');
        subset_indices = cell(1, size(unique_values, 1));
        
        % Group trials based on unique values in behavioral data
        for i = 1:size(unique_values, 1)
            subset_indices{i} = find(ismember(behav_data, unique_values(i, :), 'rows'));
        end

        % Shuffle each subset independently over the first dimension
        for i = 1:size(unique_values, 1)           
            temp_trials = shuffled_data(subset_indices{i}, :, :);
            for col = 1:size(temp_trials, 2)
                temp_trials(:, col, :) = temp_trials(randperm(length(subset_indices{i})), col, :);
            end
            shuffled_data(subset_indices{i}, :, :) = temp_trials;
        end
    end
end

% Shuffle dimensions based on index vector
for idx = 2:length(index)
    if index(idx) == 1
        if idx == 2
            warning('Shuffling the second dimension (neuralObjects) may alter probability distributions of individual objects.');
        end
        
        % Permute dimensions for shuffling
        perm_order = [idx, 1:idx-1, idx+1:nDims];
        shuffled_data_tmp = permute(shuffled_data, perm_order);
        
        % Compute size vector for the permuted data
        sizeVector = zeros(1, nDims);
        for i = 1:nDims
            sizeVector(i) = size(shuffled_data_tmp, i);
        end
        
        % Randomly permute indices along the first dimension
        [~, ind ] = sort(rand(sizeVector));
        indices = cell(1, numel(sizeVector));

        for i = 1:numel(sizeVector)
            indices{i} = 1:sizeVector(i);
        end

        [indices{:}] = ndgrid(indices{:});
        indices(1) = [];
       
        sub2ind_function = @(varargin) sub2ind(sizeVector, varargin{:});     
        shuffled_data_tmp_2 = shuffled_data_tmp(sub2ind_function(ind, indices{:}));
        
        % Invert the permutation order and update the shuffled data
        inv_perm_order = zeros(size(perm_order));
        inv_perm_order(perm_order) = 1:numel(perm_order);
        shuffled_data = permute(shuffled_data_tmp_2, inv_perm_order);
    end
end
end
