function [corrected_v, naive_v, shuff_all] = correction(inputs, outputs, corr, corefunc, varargin)
%%% *function [corrected_v, naive_v, shuff_all] = correction(inputs, outputs, corr, corefunc, opts)*
%%%
%%% The `correction` function computes bias-corrected values for the information measures implemented in the 
%%% MINT toolbox, using various bias correction methods. The function allows users to select different correction 
%%% techniques to improve the accuracy of estimated values.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input data. The data may represent
%%%             neural activity or other time-dependent signals, structured as:
%%%             - nDims [X nTimepoints] X nTrials (can be 2D or 3D).
%%%
%%%   - outputs: A corresponding matrix or array with the corrected results
%%%   
%%%
%%%   - corr: A string specifying the correction method to be used. Possible values include:
%%%           - 'qe'         : Quadratic extrapolation bias correction.
%%%           - 'le'         : Linear extrapolation bias correction.
%%%           - 'qe_shuffSub': Quadratic extrapolation combined with shuffle subtraction.
%%%           - 'le_shuffSub': Linear extrapolation combined with shuffle subtraction.
%%%           - 'shuffSub'   : Shuffle subtraction bias correction.
%%%           - 'pt'         : Panzeri-Treves bias correction.
%%%           - 'bub'        : BUB bias correction.
%%%           - 'shuffCorr'  : Shuffle correction.
%%%
%%%   - corefunc: A function call representing the core function to be used for bias correction 
%%%               (e.g., @TE for Transfer Entropy).
%%%
%%%   - varargin: Optional arguments as a structure (`opts`), including fields:
%%%              - bias: The bias correction method (same as `corr` above).
%%%              - shuff: Number of shuffling iterations (for 'shuffSub' and 'qe_shuffSub').
%%%              - xtrp: Number of extrapolations (for 'qe' and 'le').
%%%              - isBinned: Whether data has already been binned (default: false).
%%%              - suppressWarnings: Boolean flag to suppress warnings (default: false).
%%%
%%% Outputs:
%%%   - corrected_v: The bias-corrected value computed by the specified correction method.
%%%   - naive_v: The naive (uncorrected) value.
%%%   - shuff_all: The result of shuffling correction (0 if not applicable).
%%%
%%% Example:
%%%   To apply shuffle subtraction bias correction to neural data, call the function as:
%%%   opts.bias = 'shuffSub'; 
%%%   [corrected_v, naive_v, shuff_all] = correction(inputs, outputs, 'shuffSub', @TE, opts);
%%%
%%% Note:
%%% Users can add their own custom correction functions in the 'tools' folder.
%%% This function is flexible for different types of input time series data, such as neural population
%%% activities or other dynamical systems.
%%%
%%% Correction Methods:
%%% The function supports several correction methods, such as quadratic and linear extrapolation, to compute
%%% biases in the measures. Shuffle-based correction techniques involve random permutations of the
%%% input data, followed by the subtraction of the mean shuffle values from the naive results.
%%%
%%% Custom Correction Strategies:
%%% Users can define their own correction strategies by creating a MATLAB function that implements
%%% a custom correction method. The function should follow the structure of the existing core function
%%% and should accept the same input/output format. To use a custom correction, the function should be
%%% placed in the `tools` directory, and the correction function should be invoked with the name of the
%%% custom strategy.
%%%
%%% Example of Custom Correction:
%%% To create a custom correction, define a new MATLAB function:
%%%
%%%     function [corrected_v, naive_v] = custom_correction(inputs, outputs, opts)
%%%         % Implement custom correction logic here.
%%%         naive_v = feval(@my_corefunc, inputs, outputs, opts);  % Compute naive estimates.
%%%         corrected_v = your_custom_logic(naive_v);  % Apply custom correction logic.
%%%     end
%%%
%%% Save this function as `custom_correction.m` in the `tools` folder. To use it, call:
%%%
%%%     [corrected_v, naive_v] = correction(inputs, outputs, 'custom_correction', @my_corefunc, opts);
%%%
%%% This allows users to extend the correction process with tailored methods for specific analyses.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%         Step 1: Check inputs, Fill missing opts with default values           %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 5
    msg = 'not enough input arguments. See `help TE` for usage info';
    error('correction:notEnoughInput', msg);
end

opts = varargin{1};
if strcmp(opts.bias, 'shuffSub')
    default_opts.shuff = 20;
    default_opts.parallel = 0;
    default_opts.shuff_cond = false;
end
if ismember(opts.bias, {'qe', 'le'})
    default_opts.xtrp = 10;
    default_opts.parallel = 0;
end
if ismember(opts.bias, {'qe_shuffSub', 'le_shuffSub'})
    default_opts.xtrp = 10;
    default_opts.shuff = 20;
    default_opts.parallel = 0;
    default_opts.shuff_cond = false;
end

if ismember(opts.bias, {'shuffCorr', 'pt', 'bub'})
    default_opts.parallel = 0;
    default_opts.shuff_cond = false;
end


if ~isfield(opts, 'supressWarnings')
    opts.supressWarnings = false;       % Default is not supress warnings
end
if ~isfield(opts, 'isBinned')
    opts.isBinned = false;              % Default assumption data is not binned
end
default_fields= fieldnames(default_opts);
is_field_present = ismember(default_fields, fieldnames(opts));
missing_fields = default_fields(~is_field_present);
for i=1:size(missing_fields,1)
    missing_field_name = missing_fields{i};
    opts.(missing_fields{i}) = default_opts.(missing_fields{i});
    if ~opts.supressWarnings
        if iscell(default_opts.(missing_field_name))
            numericValue = cell2mat(default_opts.(missing_field_name));
        else
            numericValue = default_opts.(missing_field_name);
        end
        fprintf('Option "%s" was not specified. Using default value "%s".\n', missing_field_name, mat2str(numericValue));
    end
end

if nargin > 5
    inputs_nD = varargin{2};
end
if ~opts.isBinned
    inputs_b = binning(inputs, opts);
    opts.isBinned = true;
else
    inputs_b = inputs;
end
shuff_all = 0;
switch corr
    case {'qe','le' 'qe_shuffSub', 'le_shuffSub'}        
        if strcmp(corr, 'qe_shuff') || strcmp(corr, 'le_shuffSub')
            [corrected_v, naive_v, shuff_all] = extrapolation(inputs_b, outputs, corr, corefunc, opts);
        else
            [corrected_v, naive_v] = extrapolation(inputs_b, outputs, corr, corefunc, opts);
        end
    case 'shuffSub'
        [corrected_v, naive_v, shuff_all] = shuffle_subtraction(inputs_b, outputs, corefunc, opts);
    case 'pt'
        [corrected_v, naive_v] = pt(inputs_b, outputs, corefunc, opts);
    case 'bub'
        [corrected_v, naive_v] = bub2(inputs_b, outputs, corefunc, opts);
    case 'shuffCorr'
        [corrected_v, naive_v] = shufflecorrection(inputs_b, outputs, corefunc, opts);
    otherwise
        tools_dir = fullfile(pwd, 'src/tools');
        func_name = strcat(corefunc, '.m');
        function_path = fullfile(tools_dir, func_name);
        if exist(function_path, 'file') == 2
            % You can change the outputs and inputs_b of your function here,
            % but make sure that it fits the output definition of this
            % function.
            [corrected_v, naive_v] = feval(func_name, inputs_b, outputs, opts);
        else
            available_functions = {'qe', 'le', 'qe_shuffSub', 'le_shuffSub', 'shuffSub', 'pt', 'bub'};
            error(['The function "%s" is not defined in the tools folder.\n', ...
                'Available options are: %s.\n', ...
                'You can define additional functions by defining your own correction function in a "%s.m" file in the tools folder.'], ...
                func_name, strjoin(available_functions, ', '), func_name);
        end
end

