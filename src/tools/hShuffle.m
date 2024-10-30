function  inputs_sh = hShuffle(inputs, varargin)
% *inputs_sh = hShuffle(inputs, reqOutputs, opts)*
% SHUFFLE - Shuffle neural data with optional consistency constraints
%
% This function shuffles neural data across specified dimensions, with options
% to enforce consistency constraints based on behavioral or other conditional inputs.
% The shuffling can be performed across trials, timepoints, or objects, and the behavior
% can be modified through additional optional parameters.
%
% Inputs:
%   - inputs: A cell array containing neural data, where each element represents a variable (e.g., 'A', 'B') 
%             with dimensions such as nTrials x nObjects x nTimepoints.
%
%   - varargin: Optional input arguments, which can include:
%       - reqOutputs: A cell array of strings specifying the shuffle rules.
%                 Each string defines the variables to be shuffled and, optionally,
%                 which variables they should be conditioned on. For example:
%                 - 'AB_C': Shuffle variables A and B conditioned on C.
%                 - 'AB': Shuffle variables A and B without conditioning.
%                 Default is {'A'}, meaning only variable 'A' will be shuffled.
%
%       - opts: A structure of additional shuffle configuration options, which can include:
%               - 'dim_shuffle': A cell array specifying which dimensions to shuffle. Possible values:
%                   - 'Trials': Shuffle across trials (default).
%                   - 'Objects': Shuffle across objects.
%                   - 'Timepoints': Shuffle across timepoints.
%               - 'suppressWarnings': Boolean flag to suppress warning messages (default: false).
%
%       - cond_input: Behavioral or external data used for conditional consistency constraints (optional).
%                     This input is optional and can be used when shuffling should be conditioned on
%                     specific data that is not included in inputs. If `cond_input` is provided, it allows for conditioning on
%                     binned or external variables.
%                     For example, you might have binned behavioral data and want to shuffle neural data
%                     in `inputs` while maintaining consistency with this binned input.
%
% Output:
%   - inputs_sh: A cell array where each element contains the shuffled neural data based on the respective rule
%                in `reqOutputs`. Each element of `inputs_sh` is a modified version of the original input, shuffled
%                according to the specified rules.
%
% Description:
% This function enables flexible shuffling of neural data across various dimensions, either with or without
% consistency constraints. Shuffling can be conditioned on specific variables to maintain dependencies in the data,
% or it can be performed without conditions. The dimensions over which the shuffling occurs are customizable, allowing
% for shuffling across trials, objects, or timepoints depending on the experimental design.
%
% Examples:
% 1. To shuffle neural data 'A' and 'B', conditioned on 'C':
%      inputs_sh = hShuffle({A, B, C}, {'AB_C'}, opts);
%
% 2. To shuffle neural data across trials without any conditioning:
%      inputs_sh = hShuffle({A, B}, {'AB'}, opts);
%
% 3. To shuffle neural data 'A', maintaining consistency based on behavioral data:
%      cond_input = ... % Behavioral data
%      inputs_sh = hShuffle({A}, {'A'}, opts, cond_input);
%
% Here, the options can be defined in `opts`, and `cond_input` applies consistency constraints for shuffling.

% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.


if nargin < 1
    msg = 'not enough input arguments. See `help shuff_test` for usage info';
    error('shuffle:notEnoughInput', msg);
end
if isempty(varargin)
    reqOutputs = {'A'};
    opts = struct();
    cond_input = 0;
elseif length(varargin)==1
    if iscell(varargin{1})
        reqOutputs = varargin{1};
        opts = struct();
    else
        opts = varargin{1};
        reqOutputs = {'A'};
    end
    cond_input = 0;
elseif length(varargin)==2
     reqOutputs = varargin{1};
     opts = varargin{2};
     cond_input = 0;
else
    reqOutputs = varargin{1};
    opts = varargin{2};
    cond_input = varargin{3};
end

if ~isfield(opts, 'supressWarnings')
    opts.supressWarnings = false;
end
if ~isfield(opts, 'dim_shuffle')
    opts.dim_shuffle = {'Trials'};
end 

nVars = length(inputs);
inputs_sh = cell(1, length(reqOutputs));

for ruleIdx = 1:length(reqOutputs)
    rule = reqOutputs{ruleIdx};  % e.g., 'AB_C', 'AB', 'A_C'
    
    % Split the rule into shuffle and condition parts
    ruleParts = strsplit(rule, '_');
    shuffleVars = ruleParts{1};  % Variables to shuffle, e.g., 'AB'
    
    if length(ruleParts) > 1
        condVars = ruleParts{2};  % Condition variables, e.g., 'C'
    else
        condVars = '';  % No conditioning
    end

    % Determine which variables to shuffle
    shuffle_var_idxs = ismember(arrayfun(@(x) char('A' + (x-1)), 1:nVars, 'UniformOutput', false), cellstr(shuffleVars'));
    
    % If there's a condition variable, collect its data
    if ~isempty(condVars)
        cond_var_idxs = ismember(arrayfun(@(x) char('A' + (x-1)), 1:nVars, 'UniformOutput', false), cellstr(condVars'));
        if cond_input == 0
            cond_input_data = inputs(cond_var_idxs);  % Extract conditioned variable data
        else
            cond_input_data = cond_input(cond_var_idxs); 
        end 
    else
        cond_input_data = [];
    end

    % % Check that dimensions match across all shuffling variables
    dimsA = size(inputs{find(shuffle_var_idxs, 1)});
    % for varIdx = find(shuffle_var_idxs)
    %     if ~isequal(size(inputs{varIdx})(2:end), dimsA)
    %         error('All shuffle variables must have the same dimensions.');
    %     end
    % end
    
    nDims = length(dimsA);
    shuffle_dim = zeros(1, nDims);  

    % Set which dimensions should be shuffled based on opts.dim_shuffle
    for i = 1:length(opts.dim_shuffle)
        dim_shuffle_entry = opts.dim_shuffle{i};
        if strcmp(dim_shuffle_entry, 'Objects')
            shuffle_dim(1) = 1;
        elseif strcmp(dim_shuffle_entry, 'Timepoints')
            shuffle_dim(2) = 1;
        elseif strcmp(dim_shuffle_entry, 'Trials')
            shuffle_dim(end) = 1;
        end
    end

    % Shuffling process
    if nDims > 2
        nTimepoints = dimsA(2);
        nTrials = dimsA(3);
    else
        nTrials = dimsA(2);
        nTimepoints = 1;
    end

    shuffled_data = inputs; 
    
    if nDims == 2
        % For 2D inputs
        for varIdx = find(shuffle_var_idxs)
            if shuffle_dim(end) == 1
                if ~isempty(cond_input_data)
                    % Conditioned shuffling
                    shuffled_data{varIdx} = shuffle_core(cond_input_data{1}', inputs{varIdx}', 1, 1)';
                else
                    % Simple shuffling
                    shuffled_data{varIdx} = shuffle_core(0, inputs{varIdx}', 0, 1)';
                end
            end
        end
    elseif nDims == 3 && nTimepoints == 1
       for varIdx = find(shuffle_var_idxs)
            inputs_data = inputs{varIdx};
            inputs_data = squeeze(inputs_data);
            if size(inputs_data,1) ~= nTrials
                inputs_data = inputs_data';
            end
            if shuffle_dim(end) == 1
                if ~isempty(cond_input_data)
                    % Conditioned shuffling
                    shuffled_data{varIdx} = shuffle_core(cond_input_data{1}', inputs_data, 1, 1)';
                else
                    % Simple shuffling
                    shuffled_data{varIdx} = shuffle_core(0, inputs_data, 0, 1)';
                end
            end
       end
    else 
        % For higher dimensional inputs
        for varIdx = find(shuffle_var_idxs)
            inputs_data = inputs{varIdx};
            dims = size(inputs_data);
            shuffled_var_data = inputs{varIdx};  % Initialize the shuffled data for this variable
            if length(dims)< 3
                if ~isempty(cond_input_data)
                    shuffled_data{varIdx} = shuffle_core(cond_input_data{1}', inputs_data, 1, 1)';
                else
                    shuffled_data{varIdx} = shuffle_core(0, inputs_data, 0, 1)';
                end
                if size(shuffled_data{varIdx}, 2) ~= nTrials
                    shuffled_data{varIdx} = shuffled_data{varIdx}';  % Transponiere, wenn notwendig
                end
            else
                for tP = 1:nTimepoints
                    data_tP = squeeze(inputs_data(:, tP, :));
                    data_tP = reshape(data_tP, dims(3), dims(1));  % Reshape for shuffling

                    if ~isempty(cond_input_data)
                        % Handle conditioned shuffle
                        if length(size(cond_input_data{1})) > 2
                            cond_input_tp = squeeze(cond_input_data{1}(:, tP, :));
                            cond_input_tp = reshape(cond_input_tp, size(cond_input_data{1}', 3), size(cond_input_data{1}, 1));
                        else
                            cond_input_tp = cond_input_data{1};
                        end
                        shuffled_tmp = shuffle_core(cond_input_tp', data_tP, 1, 1)';
                    else
                        % Simple shuffle
                        shuffled_tmp = shuffle_core(0, data_tP, 0, 1)';
                    end

                    shuffled_var_data(:, tP, :) = shuffled_tmp;  % Assign shuffled data
                end
                 shuffled_data{varIdx} = shuffled_var_data;
            end          
        end
    end

    % Store the shuffled data for the current rule
    inputs_sh{ruleIdx} = shuffled_data;
end
if length(reqOutputs) == 1
    inputs_sh = inputs_sh{1}; 
end
end