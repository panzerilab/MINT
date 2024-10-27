function [ bitMask ] = clusterStatistics(infoMatrix, infoMatrixSh, significanceThreshold, clusterPercentilThreshold, pool, varargin)
% clusterStatistics - Compute cluster-based statistics for significance testing of information quantities.
% 
% This function performs a cluster-based statistical test on an input matrix of information quantities
% ('infoMatrix') by comparing it to a null distribution ('infoMatrixSh') derived from shuffled data. It uses 
% cluster-forming thresholds and a specified significance level to determine which clusters of information 
% values are statistically significant.
%
% Inputs:
%   - infoMatrix: A 2D rectangular matrix of information quantity values. This matrix represents information
%               quantities over different spatial or temporal locations.
%
%   - infoMatrixSh: A 3D matrix of shuffled information quantities, where the first two dimensions match 
%                 'infoMatrix', and the third dimension represents the number of shuffling iterations. 
%                 This is used to compute the null distribution.
%
%   - significanceThreshold: A scalar value defining the threshold for statistical significance. Clusters
%                            will be deemed significant if they survive a comparison to this threshold.
%
%   - clusterPercentilThreshold: A scalar between 0 and 1 defining the cluster-forming threshold. Only
%                                elements of 'infoMatrix' that exceed the specified percentile of the null
%                                distribution (from 'infoMatrixSh') will be considered part of a cluster.
%
%   - pool: A flag (1 or 0) indicating how the null distribution should be used:
%           - If 1, the null values are pooled across all samples in 'infoMatrix'.
%           - If 0, sample-specific null hypotheses are used, meaning each sample in 'infoMatrix' has its 
%             own corresponding null distribution.
%
%   - varargin: (Optional) Specifies the number of neighbors to consider for cluster connectivity. 
%               If not provided, the default value is 8 (8-neighbor connectivity for 2D matrices).
%
% Outputs:
%   - bitMask: A logical matrix of the same size as 'infoMatrix'. Elements with a value of 1 correspond 
%              to statistically significant clusters, while 0 indicates non-significant areas.
%
% Example:
% Suppose you have a 2D matrix of information quantities 'infoMatrix', and a corresponding 3D matrix 
% of shuffled data 'infoMatrixSh'. You want to compute a cluster-based statistical test using a 99th percentile 
% threshold for cluster formation and a significance level of 0.01, without pooling across samples:
%
%   mask = clusterStat_pool(infoMatrix, infoMatrixSh, 0.01, 0.99, 0);
%
% This will return a binary mask ('mask') where significant clusters of information quantities are marked 
% as 1, based on the given cluster-forming threshold and significance level.
%
% Additional Details:
% The function uses the `createClusters` function to identify clusters of information quantities that exceed
% the cluster-forming threshold. These clusters are compared to clusters formed in the null distribution 
% (from the shuffled data) to determine statistical significance.
% 
% Reference:
% This procedure is based on the method described in Combrisson et al.
% (2022), NeuroImage.

% This code is written by Hamed Nili (https://scholar.google.co.uk/citations?user=QBgOje0AAAAJ&hl=en)
% and Marco Celotto (https://scholar.google.it/citations?user=6nfM7V0AAAAJ&hl=it)

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


if isempty(varargin)
    numNeighbour = 8;
else 
    numNeighbour = varargin{1};
end 
threshold = significanceThreshold;

if pool % pool null values across infoMatrix dimensions (assuming that the null distribution is the same)
    allVals = infoMatrixSh(:);
    clusterThresh = prctile(allVals,clusterPercentilThreshold*100);    
    % create clusters of interest from the given 2d matrix
    [L, num, clusterVals] = createClusters(infoMatrix, clusterThresh,numNeighbour);
else % use sample-specific null hypothesis for each sample in infoMatrix (if the null distribution changes across samples)
   clusterThresh = prctile(infoMatrixSh,clusterPercentilThreshold*100,3);
        
   % create clusters of interest from the given 2d matrix
   [L, num, clusterVals] = createClusters(infoMatrix, clusterThresh,numNeighbour);
end

% count, how many times there is a cluster in the shuffled values
% with higher value than the previously computed clusters of interest
higherClusterValCount = zeros(size(clusterVals));
for s = 1:size(infoMatrixSh, 3)
    [~, ~, surrogateClusterVals] = createClusters(squeeze(infoMatrixSh(:, :, s)), clusterThresh,numNeighbour);
    if isempty(surrogateClusterVals)
        continue
    end
    higherClusterValCount(clusterVals < max(surrogateClusterVals)) = higherClusterValCount(clusterVals < max(surrogateClusterVals)) + 1;
    %maxclus(s) = max(surrogateClusterVals);
end

% put all the significant clusters into one bit mask
bitMask = zeros(size(infoMatrix));
for c = 1:num
    if (higherClusterValCount(c) / size(infoMatrixSh, 3) < threshold)
        bitMask(L == c) = 1;
    end
end

bitMask = logical(bitMask);

