function [R_1d, resps_grid] = reduce_dim(R_Nd, dim_to_collapse)
% *function [R_1d, resps_grid] = reduce_dim(R_Nd, dim_to_collapse)*
%
% Collapses the specified dimension into a single index dimension by
% enumerating ALL possible patterns (Cartesian products) formed by the
% unique values taken along the collapsed dimension. The index at each
% remaining position points to the corresponding row in resps_grid.
%
% Inputs:
%   - R_Nd : N-D array (N >= 2)
%   - dim_to_collapse : dimension to collapse (default = 1)
%
% Outputs:
%   - R_1d : array of size [1, size(R_Nd, setdiff(1:end, dim_to_collapse))]
%            containing indices into resps_grid
%   - resps_grid : M x C matrix listing ALL Cartesian patterns, where
%                  C = size(R_Nd, dim_to_collapse) and
%                  M = prod(cellfun(@numel, resps)) is the number of combos
%
% Note:
%   Every observed pattern in R_Nd must appear in resps_grid because
%   resps_grid is built from the per-component uniques (Cartesian product).
%
% EXAMPLE:
%   R_Nd = randi(3,[4,5,3]);                % values 1..3
%   [R_1d, resps_grid] = reduce_dim(R_Nd,1); % collapse dim 1 (size 4)
%
% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT. Licensed under GPLv3 or later.

    if nargin < 1
        error('reduce_dim:notEnoughInput','Please input your data.');
    end
    if nargin < 2
        warning('reduce_dim:notEnoughInput', ...
            'Not enough input arguments. Collapsing the first dimension by default.');
        dim_to_collapse = 1;
    end
    if any(isnan(R_Nd(:)))
        error('reduce_dim:NaNInput','R_Nd contains NaNs. Aborting.');
    end

    dims = size(R_Nd);
    Ndims = numel(dims);

    % Trivial case: nothing to collapse
    if dims(dim_to_collapse) <= 1
        R_1d = R_Nd;
        resps_grid = R_Nd(:).'; % degenerate listing
        return;
    end

    % Bring the collapsed dimension to the front for easy reshaping
    remaining_dims = setdiff(1:Ndims, dim_to_collapse, 'stable');
    R_perm = permute(R_Nd, [dim_to_collapse, remaining_dims]);

    C = size(R_perm,1);                % length of collapsed dimension
    rest_shape = size(R_perm);         % [C, d2, d3, ...]
    rest_shape = rest_shape(2:end);    % [d2, d3, ...]
    K = prod(rest_shape);              % number of positions in remaining dims

    % Flatten so each column corresponds to one position in the remaining dims
    % Patterns are rows of size C (one value per collapsed component)
    R2 = reshape(R_perm, C, K);   % C x K
    patterns = R2.';              % K x C

    % Unique values per component (across all positions), then Cartesian product
    resps = cell(1, C);
    for c = 1:C
        resps{c} = unique(R2(c, :));   % sorted uniques per component
    end

    % Build Cartesian product grid: M x C
    [grid_cells{1:C}] = ndgrid(resps{:});
    resps_grid = reshape(cat(C+1, grid_cells{:}), [], C);  % M x C

    % Map every observed pattern to its row in the Cartesian grid
    % Since resps_grid was built from per-component uniques, every pattern must match.
    [~, idx] = ismember(patterns, resps_grid, 'rows');     % K x 1

    % Reshape to output size [1, rest_shape]
    R_1d = reshape(idx, [1, rest_shape]);

end
