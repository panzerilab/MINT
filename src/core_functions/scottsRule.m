function num_bins = scottsRule(data)   
% Description:
% Calculates the number of bins using Scott's rule for histogram bin width.
%
% Inputs:
% - data: Input data vector.
%
% Outputs:
% - num_bins: Optimal number of bins calculated using Scott's rule.
%
% Source:
% This implementation is based on Scott's rule for histogram bin width.
% Reference: Scott, D. W. (1979). On optimal and data-based histograms. Biometrika, 66(3), 605-610.
% DOI: https://doi.org/10.1093/biomet/66.3.605

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

    N = numel(data);
    h = 3.49 * std(data) * N^(-1/3);
    data_range = range(data);
    if data_range == 0 || h == 0
        num_bins = 1;
    else
        num_bins = ceil(data_range / h);
    end

end
