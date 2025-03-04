function inputs = nan_method(inputs,handling)
% *function inputs = nan_method(inputs, handling)*
%
% The `nan_method` function manages NaN values in input data based on the specified handling method.
% This function supports three methods for handling NaN values: removing trials with NaN values,
% setting NaNs to zero, or throwing an error if NaN values are present.
%
% Inputs:
%   - inputs:   A cell array containing the input data sets. Each element represents a different
%               variable or dataset, typically organized as matrices where the dimensions represent
%               different trials and data points.
%
%   - handling: A string specifying how to handle NaN values in `inputs`. Available options are:
%               - 'removeTrial' : Removes all trials (across all variables) that contain NaN values.
%               - 'error'       : (default) Throws an error if NaN values are detected in any input.
%
% Outputs:
%   - inputs:   A cell array with NaN values processed according to the specified handling method.

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

if strcmp(handling, 'removeTrial')
    NaN_indices = [];
    for var = 1:length(inputs)
        var_data = inputs{var};
        nanIndices = any(isnan(var_data), 1:ndims(var_data) - 1);
        nanIndices_Trials = find(nanIndices);
        NaN_indices = [NaN_indices, nanIndices_Trials'];
    end
    for var = 1:length(inputs)
        var_data = inputs{var};
        dims =ndims(var_data);
        if dims == 3
            inputs{var}(:, :, NaN_indices) = [];
        elseif dims == 2
            inputs{var}(:, NaN_indices) = [];
        end
    end
% elseif strcmp(handling, 'setToZero')
%     for var = 1:length(inputs)
%         inputs{var}(isnan(inputs{var})) = 0;
%     end
elseif strcmp(handling, 'error')
    for var = 1:length(inputs)
        if any(isnan(inputs{var}),'all')
            letter = char(64 + var);
            msg = sprintf('NaN value found in input variable %s. Use opts.NaN_handling to handle NaN values (see help nan_method for more information)', letter);
            error('checkInputs:NaNDetected', msg);  
        end
    end
end
end

