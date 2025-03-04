function data_out = create_NullDistribution_groupLevel(data, null_samples)
% create_NullDistribution_groupLevel - Generate a null distribution for group-level data using bootstrapped shuffling.
%
% This function performs a bootstrap-based shuffling procedure at the group level to create a null distribution 
% of the input data. It is designed to handle multidimensional data arrays, where the first dimension represents 
% repetitions or subjects, and the last dimension corresponds to shuffles. The function reshuffles the data across 
% repetitions, computes the mean for each shuffle, and outputs a reshuffled distribution for statistical comparison.
%
% Inputs:
%   - data: A multidimensional array representing the data to analyze. The array should have at least three 
%           dimensions, where the first dimension corresponds to repetitions or subjects, and the last dimension 
%           represents shuffled samples. For example, a 3D array could be structured as (nReps, nDims, nShuffles).
%
%   - null_samples: An integer specifying the number of shuffled (null) samples to generate. This determines how 
%                   many times the data will be reshuffled to create the null distribution.
%
% Outputs:
%   - data_out: A multidimensional array containing the shuffled data averages, with the same dimensions as the 
%               input data, except for the first dimension (repetitions), which is removed. The output dimensions 
%               depend on the input data structure, and the last dimension of 'data_out' will correspond to the 
%               number of null samples generated.
%
% Example:
% Suppose we have a 3D dataset 'data', where the first dimension corresponds to 10 subjects, the second dimension 
% represents 50 features, and the third dimension contains 100 shuffled samples. To create a null distribution 
% with 1000 bootstrapped samples:
%
%   nullDist = create_NullDistribution_groupLevel(data, 1000);
%
% The resulting 'nullDist' will have dimensions (50, 1000), where 50 corresponds to the original features, and 1000 
% represents the null samples generated through bootstrapping.
%
% Additional Details:
% The function adapts to the dimensionality of the input data (3D to 6D). It reshuffles the last dimension across 
% repetitions and computes the mean for each shuffle, creating a null distribution for statistical inference.
%
% Reference:
% This method is developed by Hamed Nili (https://scholar.google.co.uk/citations?user=QBgOje0AAAAJ&hl=en)

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

nReps = size(data,1);
nDim = numel(size(data));

size_surrogate = size(data);
nShuff = size_surrogate(end);
size_surrogate(end) = [];

size_out = size(data);
size_out(1) = [];

data_out = nan(size_out);

for bIdx = 1:null_samples
    randIdxs = randi(nShuff,1,nReps); 
    tmpSurrogate = nan(size_surrogate);
    if nDim == 3
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:) = data(repIdx,:,randIdxs(repIdx));
        end
        data_out(:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 4
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:) = data(repIdx,:,:,randIdxs(repIdx));
        end
        data_out(:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 5
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:,:) = data(repIdx,:,:,:,randIdxs(repIdx));
        end
        data_out(:,:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    elseif nDim == 6
        for repIdx = 1:nReps
            tmpSurrogate(repIdx,:,:,:,:) = data(repIdx,:,:,:,:,randIdxs(repIdx));
        end
        data_out(:,:,:,:,bIdx) = squeeze(mean(tmpSurrogate,1));
    end
end

end

