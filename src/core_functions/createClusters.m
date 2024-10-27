function [ L, num, values ] = createClusters(infQuant, clusterThreshold,numNeighbour)
% createClusters - Identify and create clusters of information quantities based on a threshold.
%
% This function processes a matrix of information quantities, identifying clusters of values
% that are above a specified threshold. It uses connected component labeling to form clusters
% based on spatial proximity (neighbor connectivity) and computes the sum of values for each cluster.
%
% Inputs:
%   - infQuant: A rectangular matrix of information quantity values. This matrix represents
%               information quantities across different spatial or temporal locations.
%
%   - clusterThreshold: A scalar value that defines the threshold for cluster formation.
%                       Only elements in the matrix greater than or equal to this threshold
%                       will be considered part of a cluster.
%
%   - numNeighbour: A scalar that defines the connectivity used for clustering.
%                   Typically, it can be 4 or 8 for 2D matrices, specifying the number of
%                   neighboring elements required to consider a point part of a cluster.
%                   For 3D matrices, this number could be higher, like 6 or 26, depending
%                   on the desired connectivity.
%
% Outputs:
%   - L: A matrix of the same size as 'infQuant', where each element is labeled by its cluster
%        number. Elements not part of any cluster are labeled as 0.
%
%   - num: The total number of clusters identified in 'infQuant' based on the threshold and 
%          connectivity.
%
%   - values: A vector of length 'num', where each entry corresponds to the sum of the 
%             information quantity values for the corresponding cluster in 'infQuant'.
%
%
% Additional Info:
% infQuant should be a rectangular matrix of values of the information
% quantity
%
% take only non-zero values beacuse sometimes we don't compute the whole
% matrix and our quantities shouldn't be exactly zero even if they are very
% low

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

allValues = infQuant(:);

% just a corner case when all values are smaller than the clusterThreshold
if max(allValues) < clusterThreshold
    L = zeros(size(infQuant));
    num = 0;
    values = [];
else
    % create a bit mask of all those that are above the threshold
    minValid = clusterThreshold;
    valid = infQuant >= minValid;
    
    % create clusters
    [L, num] = bwlabeln(valid,numNeighbour); % you can set the second parameter to 4 if you want clusters to be connected only with the 4 closest neighbours, or 8 to take all around

    % sum the t-values for them
    values = zeros(num, 1);
    for i = 1:num
        values(i) = sum(sum(infQuant(L == i)));
    end
end

