function [corrected_v, naive_v, shuff_all] = correction(inputs, reqOutputs, corr, corefunc, varargin)
% correction - Compute bias-corrected and naive information-theoretic values
%
% The `correction` function computes bias-corrected values for the information measures implemented in the 
% MINT toolbox, using various bias correction methods. The function allows users to select different correction 
% techniques to improve the accuracy of estimated values.
%
% Inputs:
%   - inputs: A cell array containing the input data. The data may represent
%             neural activity or other time-dependent signals, structured as:
%             - nDims [X nTimepoints] X nTrials (can be 2D or 3D).
%
%   -  reqOutputs: A cell array of strings specifying which measures to compute
%   
%
%   - corr: A string specifying the correction method to be used. Possible values include:
%           - 'qe'         : Quadratic extrapolation bias correction.
%           - 'le'         : Linear extrapolation bias correction.
%           - 'qe_shuffSub': Quadratic extrapolation combined with shuffle subtraction.
%           - 'le_shuffSub': Linear extrapolation combined with shuffle subtraction.
%           - 'shuffSub'   : Shuffle subtraction bias correction.
%           - 'pt'         : Panzeri-Treves bias correction.
%           - 'bub'        : best upper bound correction.
%           - 'shuffCorr'  : Shuffle correction.
%
%   - corefunc: A function call representing the core function to be used for bias correction 
%               (e.g., @TE for Transfer Entropy).
%
%   - varargin: Optional arguments as a structure (`opts`), including fields:
%              - bias: The bias correction method (same as `corr` above).
%              - shuff: Number of shuffling iterations (for 'shuffSub' and 'qe_shuffSub').
%              - xtrp: Number of extrapolations (for 'qe' and 'le').
%              - isBinned: Whether data has already been binned (default: false).
%              - suppressWarnings: Boolean flag to suppress warnings (default: false).
%
% Outputs:
%   - corrected_v: The bias-corrected value computed by the specified correction method.
%   - naive_v: The naive (uncorrected) value.
%   - shuff_all: The result of shuffling correction (0 if not applicable).
%
% Example:
%   To apply shuffle subtraction bias correction, call the function as:
%   opts.bias = 'shuffSub'; 
%   [corrected_v, naive_v, shuff_all] = correction(inputs, reqOutputs, 'shuffSub', @TE, opts);
%
% Note:
% Users can add their own custom correction functions in the 'tools' folder.
% This function is flexible for different types of input time series data, such as neural population
% activities or other dynamical systems.
%
% Correction Methods:
% The function supports several correction methods, such as quadratic and linear extrapolation, to compute
% biases in the measures. Shuffle-based correction techniques involve random permutations of the
% input data, followed by the subtraction of the mean shuffle values from the naive results.
%
% Custom Correction Strategies:
% Users can define their own correction strategies by creating a MATLAB function that implements
% a custom correction method. The function should follow the structure of the existing core function
% and should accept the same input/output format. To use a custom correction, the function should be
% placed in the `tools` directory, be added to the path, and the correction function should be invoked with the name of the
% custom strategy.
%
% Save this function as `custom_correction.m`. To use it, call:
%     [corrected_v, naive_v] = correction(inputs, reqOutputs, 'custom_correction', @my_corefunc, opts);
% This allows users to extend the correction process with tailored methods for specific analyses.

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
%         Step 1: Check inputs, Fill missing opts with default values           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    msg = 'not enough input arguments. See `help correction` for usage info';
    error('correction:notEnoughInput', msg);
end

opts = varargin{1};
if strcmp(corr, 'shuffSub')
    defaultOpts.shuff = 20;
    defaultOpts.parallel = 0;
elseif ismember(corr, {'qe', 'le'})
    defaultOpts.xtrp = 10;
    defaultOpts.parallel = 0;
elseif ismember(corr, {'qe_shuffSub', 'le_shuffSub'})
    defaultOpts.xtrp = 10;
    defaultOpts.shuff = 20;
    defaultOpts.parallel = 0;
elseif ismember(corr, {'shuffCorr', 'pt', 'bub'})
    defaultOpts.parallel = 0;
else
    defaultOpts = struct();
end

if ~isfield(opts, 'supressWarnings')
    opts.supressWarnings = false;       % Default is not supress warnings
end
if ~isfield(opts, 'isBinned')
    opts.isBinned = false;              % Default assumption data is not binned
end
default_fields= fieldnames(defaultOpts);
is_field_present = ismember(default_fields, fieldnames(opts));
missing_fields = default_fields(~is_field_present);
for i=1:size(missing_fields,1)
    missing_field_name = missing_fields{i};
    opts.(missing_fields{i}) = defaultOpts.(missing_fields{i});
    if ~opts.supressWarnings
        if iscell(defaultOpts.(missing_field_name))
            numericValue = cell2mat(defaultOpts.(missing_field_name));
        else
            numericValue = defaultOpts.(missing_field_name);
        end
        fprintf('Option "%s" was not specified. Using default value "%s".\n', missing_field_name, mat2str(numericValue));
    end
end

if nargin > 5
    inputs_nD = varargin{2};
end


shuff_all = 0;

switch corr
    case {'qe','le' }
        [corrected_v, naive_v] = extrapolation(inputs, reqOutputs, corr, corefunc, opts);
    case {'qe_shuffSub', 'le_shuffSub'}
        [corrected_v, naive_v] = shuffSub_extrapolation(inputs, reqOutputs, corr, corefunc, opts);
    case 'shuffSub'
        [corrected_v, naive_v, shuff_all] = shuffle_subtraction(inputs, reqOutputs, corefunc, opts);
    case 'shuffCorr'
        [corrected_v, naive_v] = shuffle_correction(inputs, reqOutputs, corefunc, opts);
    otherwise
        func_handle = str2func(corr);
        if exist(corr, 'file') == 2
            % You can change the reqOutputs and inputs_b of your function here,
            % but make sure that it fits the output definition of this
            % function.
            [corrected_v, naive_v] = feval(func_handle, inputs, reqOutputs, opts);
        else
            available_functions = {'qe', 'le', 'qe_shuffSub', 'le_shuffSub', 'shuffSub', 'pt', 'bub'};
            msg = sprintf(['The function "%s" is not defined in the tools folder.\n', ...
                'Available options are: %s.\n', ...
                'You can define additional functions by defining your own correction function in a "%s" file in the tools folder.'], ...
                corr, strjoin(available_functions, ', '), corr);
            error('Correction:UndefinedFunction', msg);
        end
end

