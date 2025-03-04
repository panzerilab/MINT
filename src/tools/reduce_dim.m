function R_1d = reduce_dim(R_Nd, dim_to_collapse)
% *function R_1d = reduce_dim(R_Nd, dim_to_collapse)*
%
% The reduce_dim function takes a multi-dimensional array and collapses the specified 
% dimension into a single dimension, effectively reducing the overall dimensionality 
% of the array. This is particularly useful for simplifying data representation 
% while preserving unique values from the collapsed dimension.
%
% Inputs:
%   - R_Nd: A multi-dimensional array with dimensions greater than or equal to 2.
%            The array can represent various data structures, such as trials, 
%            timepoints, or other features.
%
%   - dim_to_collapse: (Optional) An integer specifying which dimension to collapse.
%                      If not provided, the first dimension is collapsed by default.
%
% Outputs:
%   - R_1d: A reduced-dimensional array where the specified dimension has been collapsed. 
%            The resulting structure will have one less dimension than the original 
%            input, maintaining the unique values from the collapsed dimension.
%
% Note: 
% The resulting array R_1d contains the indices of the unique values from the 
% collapsed dimension. This allows for efficient representation of the data while 
% enabling further analysis without losing critical information.
%
% EXAMPLE
% Suppose we have a 3D array representing responses from multiple trials, 
% with the dimensions representing trials, features, and observations:
% 
% R_Nd = randn(4, 5, 3);  
% To reduce the dimensionality by collapsing the first dimension, the function can be called as:
% R_1d = reduce_dim(R_Nd, 1);
% R_1d will then be a 1 x 5 x 3 double.
%
% Alternatively, to collapse the second dimension, you would call:
% R_1d = reduce_dim(R_Nd, 2);
% R_1d will then be a 4 x 1 x 3 double.

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
if nargin < 1
    msg = 'Please input your data.';
    error('reduce_dim:notEnoughInput', msg);
end

if nargin < 2
    warning('reduce_dim:notEnoughInput', 'Not enough input arguments. Collapsing the first dimension by default.');
    dim_to_collapse = 1;
end
if any(any(isnan(R_Nd)))
    msg = "R_Nd contains NaNs. Aborting.";
    error('reduce_dim:NaNInput', msg);
end
if size(R_Nd, dim_to_collapse) > 1
    dims = size(R_Nd);
    Ndims = length(dims);
    remaining_dims = setdiff(1:length(dims), dim_to_collapse);
    R_Nd = permute(R_Nd, [dim_to_collapse, remaining_dims]);
    dims = size(R_Nd);
    collapsed_dim_size = size(R_Nd, 1);
    resps = cell(1, collapsed_dim_size);
    for d = 1:collapsed_dim_size
        resps{d} = unique(R_Nd(d,:));
    end
    [resps_grid{1:collapsed_dim_size}] = ndgrid(resps{:});
    resps_grid = reshape(cat(collapsed_dim_size+1, resps_grid{:}), [], collapsed_dim_size);
    R_1d = zeros([1,dims(2:end)]);
    nRemainingDims = length(remaining_dims);
    indices = cell(1, nRemainingDims);
    for d = 1:nRemainingDims
        indices{d} = 1:size(R_Nd, remaining_dims(d));
    end
    [indices_grid{1:nRemainingDims}] = ndgrid(indices{:});
    indices_grid = reshape(cat(nRemainingDims+1, indices_grid{:}), [], nRemainingDims);
    if Ndims == 2
        [~, idx] = ismember(R_Nd', resps_grid, 'rows');
        R_1d(idx > 0) = idx(idx > 0);
    elseif Ndims == 3
        for k = 1:dims(2)
            slice = squeeze(R_Nd(:,k,:))';
            [~, idx] = ismember(slice, resps_grid, 'rows');
            R_1d(1,k,:) = idx;
        end
    else
        for k = 1:size(indices_grid,1)
            idx = num2cell(indices_grid(k,:));
            R_Nd_slice = R_Nd(:, idx{:});
            logical_idx = all(R_Nd_slice' == resps_grid, 2);
            loc = find(logical_idx);
            R_1d(:, idx{:}) = loc;
        end
    end
else
    R_1d = R_Nd;
end

end

