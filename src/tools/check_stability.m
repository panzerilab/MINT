function [results_full, results_partition21, results_partition22] = check_stability(inputs,corefunc,varargin)
% *function [results_full, results_partition21, results_partition22] = check_stability(inputs, corefunc, varargin)*
%
% The `check_stability` function evaluates the stability of core computations in the MINT toolbox by
% calculating specified core function (e.g., PID or other) outputs on a full dataset and two equal
% partitions of it. This split-sample approach helps in assessing the reliability and consistency of
% information measures across dataset splits.
%
% Inputs:
%   - inputs: A cell array containing data for multiple variables or sources. Each entry represents
%             a data variable (e.g., neural recordings from different regions or channels):
%             - inputs{1}, inputs{2}, ... inputs{N} : Each input variable in matrix form,
%               with dimensions nDims x (nTimepoints x) nTrials. 
%
%   - corefunc: A function handle pointing to a specific core function from the MINT toolbox (e.g., PID),
%               used to calculate information-theoretic values on the input data. For more details on 
%               specific options available for each core function, refer to 'help 'corefunc'' e.g 'help PID'.
%               MINT options:  @H, @MI, @cMI, @TE, @cTE, @PID, @II, @FIT, @cFIT, 
%               - you can also use your use your own func
%
%   - varargin: Optional arguments that are passed as a structure, typically including fields like:
%               - `bin_method`      : Binning method to apply (e.g., 'none', 'eqpop', 'eqspace', etc.)
%               - `n_bins`          : Specifies the number of bins for binning (default is {3})
%               - `bias`            : Indicates the bias correction method ('naive', 'qe', 'ShuffSub', etc.)
%               - Additional options may vary according to the core function specified
%               Note: by default this function does not compute the nulldistribution so computeNulldist is always set to false
%
% Outputs:
%   - results_full: Struct containing the output values of `corefunc` applied to the full dataset:
%       - results_full.corrected : Bias-corrected results from `corefunc`.
%       - results_full.naive     : Naive results from `corefunc` without bias correction.
%
%   - results_partition21: Struct containing the output values of `corefunc` applied to the first data partition:
%       - results_partition21.corrected : Corrected values from `corefunc` on the first partition.
%       - results_partition21.naive     : Naive values from `corefunc` on the first partition.
%
%   - results_partition22: Struct containing the output values of `corefunc` applied to the second data partition:
%       - results_partition22.corrected : Corrected values from `corefunc` on the second partition.
%       - results_partition22.naive     : Naive values from `corefunc` on the second partition.
%
% Example:
%   To evaluate stability in partial information decomposition results between two neural sources about a target,
%   the function can be called as follows:
%   inputs = {source1_data, source2_data, target_data};
%   outputs = {'PID_atoms'};
%   results = check_stability(inputs, @PID, outputs, opts);
% 
% Here, `opts` represents additional options tailored for the `corefunc` (e.g., PID options),
% and '@PID' represents any core MINT function handle selected for stability analysis.
%
% For further details on available options, consult the documentation specific to the chosen core function.

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
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses>


if nargin < 2
    msg = 'not enough input arguments. See `help check_stability` for usage info';
    error('checkStability:notEnoughInput', msg);
end
[outputs, opts] = check_inputs(corefunc,inputs,varargin);

DimsA = size(inputs{1});
nTrials = nDimsA(end);

opts.computeNulldist = false;

if ~opts.isBinned
    inputs_b = binning(inputs,opts);
    opts.isBinned = true;
end
inputs_1d = inputs_b;
for var = 1:length(inputs_b)
    sizeVar = size(inputs_1d{var});
    if sizeVar(1) > 1
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
    end
end
ri = randperm(nTrials, nTrials);
for var = length(inputs_1d)
    inputs_1d{var} = inputs_1d{var}(ri);
end

inputs_p_1 = partition(inputs_1d, 2, 1 ,1);
inputs_p_2 = partition(inputs_1d, 2, 2, 1);


[value_corr, value_naive] = feval(corefunc, inputs_1d, outputs, opts);
[value_corr_1, value_naive_1] = feval(corefunc, inputs_p_1, outputs, opts);
[value_corr_2, value_naive_2] = feval(corefunc, inputs_p_2, outputs, opts);

results_full.corrected = value_corr;
results_full.naive = value_naive;
results_partition21.corrected = value_corr_1;
results_partition21.naive = value_naive_1;
results_partition22.corrected = value_corr_2;
results_partition22.naive = value_naive_2;
end
