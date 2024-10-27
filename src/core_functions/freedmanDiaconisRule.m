function num_bins = freedmanDiaconisRule(data)
% Description:
% Calculates the optimal number of bins using the Freedman-Diaconis rule for histogram bin width.
%
% Inputs:
% - *data*: Input data vector.
%
% Outputs:
% - *num_bins*: Optimal number of bins calculated using the Freedman-Diaconis rule.
%
% Source:
% This implementation is based on the Freedman-Diaconis rule for histogram bin width.
% Reference: Freedman, D. and Diaconis, P. (1981). On the histogram as a density estimator: L2 theory.
% DOI: https://doi.org/10.1007/BF01025868

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
IQR = iqr(data); % Interquartile range
h = 2 * IQR * N^(-1/3); % Bin width based on Freedman-Diaconis rule
data_range = range(data);

if data_range == 0 || h == 0
    num_bins = 1;
else
    num_bins = ceil(data_range / h);
end

end