function null_distribution = create_nullDist(inputs, outputs, corefunc, varargin)
% *function null_distribution = create_NullDistribution(inputs, outputs, corefunc, varargin)*
%
% The `create_NullDistribution` function generates a null distribution through random shuffling
% and function evaluation for the purpose of statistical testing. This method is repeatedly shuffling
% inputs and calculating the output of the specified core function, `corefunc`.
%
% Inputs:
%   - inputs:  A cell array containing input datasets. Each element corresponds to a different
%              variable or time series. The dimensions of the inputs are context-specific and can
%              vary based on the application.
%
%   - outputs: A list of desired outputs corresponding to different evaluations or statistics computed
%              by the core function, such as entropy, transfer entropy, or other metrics.
%
%   - corefunc: A function handle that computes a specific statistic or metric from the provided
%               inputs and outputs. The core function must accept shuffled inputs and apply the
%               desired computation for generating the null distribution.
%
%   - varargin: Additional options provided as a structure. The options allow for flexibility
%               in terms of binning methods, shuffle procedures, number of samples, and parallelization.
%               Some key optional fields include:
%               - `bias`: A string specifying bias correction methods (default: 'naive').
%               - `bin_method`: A cell array for binning methods (default: {'none'}).
%               - `n_bins`: The number of bins for discretizing continuous data (default: 3).
%               - `shuff`: An integer specifying the number of shuffles to apply (default: 0).
%               - `n_samples`: Number of null samples to generate (default: 100).
%               - 'shuffling': Cell with string specifiying what to shuffle (can be also more than one)(for help type "help hShuffle" function)(default: {'A'})
%               - `parallel_sampling`: Boolean for enabling parallel sampling (default: false).
%               - `dim_shuffle`: Specifies the dimension along which shuffling occurs (default: {'Trials'}).
%               - `suppressWarnings`: Boolean flag to suppress warnings during execution (default: false).
%
% Outputs:
%   - null_distribution: A cell array containing the generated null distribution of outputs.
%                        If conditional shuffling is enabled, this output may also contain
%                        a second null distribution for unconditioned shuffling.
%
% Notes:
% This function is typically used for statistical hypothesis testing, allowing users to
% generate a distribution of random outputs to compare against the actual observed data.

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check Inputs, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_opts.bias = 'naive';
default_opts.bin_method = {'none'};
default_opts.n_bins = {3};
default_opts.n_samples = 100;
default_opts.shuffling = {'A'};
default_opts.parallel = 0;
default_opts.dim_shuffle = {'Trials'};

if nargin < 1
    msg = 'not enough input arguments. See `help create_NullDistribution` for usage info';
    error('createNullDistribution:notEnoughInput', msg);
end
opts = varargin{1};
if ~isfield(opts, 'suppressWarnings')
    opts.suppressWarnings = false;
end
default_fields= fieldnames(default_opts);
is_field_present = ismember(default_fields, fieldnames(opts));
missing_fields = default_fields(~is_field_present);
for i=1:size(missing_fields,1)
    missing_field_name = missing_fields{i};
    opts.(missing_fields{i}) = default_opts.(missing_fields{i});
    if ~opts.suppressWarnings
        if iscell(default_opts.(missing_field_name))

            numericValue = cell2mat(default_opts.(missing_field_name));
        else
            numericValue = default_opts.(missing_field_name);
        end
        fprintf('Option "%s" was not specified. Using default value "%s".\n', missing_field_name, mat2str(numericValue));
    end
end

nTimepoints_all = 1;
Dims_tmp = size(inputs{1});
nTrials_comp = Dims_tmp(end);
for var = 1:length(inputs)
    Dims_var = size(inputs{var});
    if length(Dims_var) == 2
        nTrials = Dims_var(end);
        if nTrials ~= nTrials_comp
            error('createNulldist: Inconsistent number of trials in input');
        end
    elseif length(Dims_var) == 3
        nTrials = Dims_var(end);
        if nTrials ~= nTrials_comp
            error('createNulldist: Inconsistent number of trials in input');
        end
        nTimepoints_all = [nTimepoints_all, Dims_var(2)];
    else
        error('createNulldist: Invalid input size');
    end
end
nTimepoints = max(nTimepoints_all);
opts_funcCall = opts;
opts_funcCall.compute_nulldist = false;

if length(opts.shuffling) >1
    null_distribution = struct();
end

if ~isfield(opts, 'parallel_sampling')
    opts.parallel_sampling = false;
end

for shIdx = 1:length(opts.shuffling)
    field_name = opts.shuffling{shIdx};
    null_distribution_tmp = cell(opts.n_samples, length(outputs));
    if ~opts.parallel_sampling
        for shuffIdx = 1:opts.n_samples
            shuffOutputs = {opts.shuffling{shIdx}};
            inputs_sh = hShuffle(inputs, shuffOutputs,opts);
            values = feval(corefunc, inputs_sh, outputs, opts_funcCall);
            for outIdx = 1:length(outputs)
                null_distribution_tmp{shuffIdx, outIdx} = values{outIdx};
            end
        end
    else
        parfor shuffIdx = 1:opts.n_samples
            shuffOutputs = {opts.shuffling{shIdx}};
            inputs_sh = hShuffle(inputs, shuffOutputs, opts);           
            values = feval(corefunc, inputs_sh, outputs, opts_funcCall);
            null_distribution_tmp(shuffIdx, :) = values;
        end
    end
    if nTimepoints == 1 || isequal(corefunc, @FIT) || isequal(corefunc, @cFIT) || isequal(corefunc, @TE) || isequal(corefunc, @cTE)
        null_distribution_tmp = cell2mat(null_distribution_tmp);
    end
    if length(opts.shuffling) > 1
        null_distribution.(field_name) = null_distribution_tmp;
    else
        null_distribution = null_distribution_tmp;
    end
end
end



