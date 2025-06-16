function result = mean_SEM(input, reqOutputs)
% mean_SEM computes the mean and SEM for each cell in dataCells based on reqOutputs.
% dataCells   - 1xN cell array, where each cell contains a data matrix.
% reqOutputs  - Cell array of strings specifying requested outputs (e.g., 'MeanAll', 'SEMRow').

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
% along with this program.  If not, see <http://www.gnu.org/licenses>

% List of possible output requests
possibleOutputs = {'MeanAll', 'MeanRow', 'MeanCol', 'SEMAll', 'SEMRow', 'SEMCol'};
% Validate reqOutputs against possibleOutputs
[isMember, ~] = ismember(reqOutputs, possibleOutputs);
if any(~isMember)
    nonMembers = reqOutputs(~isMember);
    msg = sprintf('Invalid reqOutputs: %s', strjoin(nonMembers, ', '));
    error('mean_SEM:invalidOutput', msg);
end

% Initialize resultCells to store outputs for each requested calculation
result = cell(size(reqOutputs));

% Process each requested output type
for idx = 1:numel(reqOutputs)
    outputType = reqOutputs{idx};
    outputData = cell(size(input));
    for i = 1:numel(input)
        data = input{i};

        switch outputType
            case 'MeanAll'  % Mean of all elements
                outputData{i} = mean(data(:));

            case 'MeanRow'  % Mean across rows (returns a column vector)
                outputData{i} = mean(data, 2);

            case 'MeanCol'  % Mean across columns (returns a row vector)
                outputData{i} = mean(data, 1);

            case 'SEMAll'  % SEM of all elements
                outputData{i} = std(data(:)) / sqrt(numel(data));

            case 'SEMRow'  % SEM across rows (returns a column vector)
                outputData{i} = std(data, 0, 2) ./ sqrt(size(data, 2));

            case 'SEMCol'  % SEM across columns (returns a row vector)
                outputData{i} = std(data, 0, 1) ./ sqrt(size(data, 1));
        end
    end
    result{idx} = outputData;
end
end
