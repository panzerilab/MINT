function null_distribution = create_NullDistribution(inputs, outputs, corefunc, varargin)
%%% *function null_distribution = create_NullDistribution(inputs, outputs, corefunc, varargin)*
%%%
%%% The `create_NullDistribution` function generates a null distribution through random shuffling
%%% and function evaluation for the purpose of statistical testing. This method is repeatedly shuffling 
%%% inputs and calculating the output of the specified core function, `corefunc`.
%%%
%%% Inputs:
%%%   - inputs:  A cell array containing input datasets. Each element corresponds to a different
%%%              variable or time series. The dimensions of the inputs are context-specific and can
%%%              vary based on the application.
%%%
%%%   - outputs: A list of desired outputs corresponding to different evaluations or statistics computed
%%%              by the core function, such as entropy, transfer entropy, or other metrics.
%%%
%%%   - corefunc: A function handle that computes a specific statistic or metric from the provided
%%%               inputs and outputs. The core function must accept shuffled inputs and apply the
%%%               desired computation for generating the null distribution.
%%%
%%%   - varargin: Additional options provided as a structure. The options allow for flexibility
%%%               in terms of binning methods, shuffle procedures, number of samples, and parallelization.
%%%               Some key optional fields include:
%%%               - `bias`: A string specifying bias correction methods (default: 'naive').
%%%               - `bin_method`: A cell array for binning methods (default: {'none'}).
%%%               - `n_bins`: The number of bins for discretizing continuous data (default: 3).
%%%               - `shuff`: An integer specifying the number of shuffles to apply (default: 0).
%%%               - `n_samples`: Number of null samples to generate (default: 100).
%%%               - `cond_var`: Conditional variables to hold constant during shuffling (default: {'C'}).
%%%               - `cond_shuff`: Boolean to control whether to perform conditional shuffling (default: false).
%%%               - `shuff_var`: Variables to shuffle (default: {'A', 'B'}).
%%%               - `parallel`: Boolean for enabling parallel processing (default: false).
%%%               - `dim_shuffle`: Specifies the dimension along which shuffling occurs (default: {'Trials'}).
%%%               - `supressWarnings`: Boolean flag to suppress warnings during execution (default: false).
%%%
%%% Outputs:
%%%   - null_distribution: A cell array containing the generated null distribution of outputs.
%%%                        If conditional shuffling is enabled, this output may also contain
%%%                        a second null distribution for unconditioned shuffling.
%%%
%%% Procedure:
%%%   1. **Input Checking and Option Handling**: The function first checks for missing optional
%%%      parameters and fills them in with default values. Warnings are shown if any options
%%%      are missing, unless `suppressWarnings` is set to true.
%%%
%%%   2. **Shuffling Process**: The inputs are shuffled according to the specified options.
%%%      The function supports both unconditional shuffling and conditional shuffling,
%%%      depending on the provided `cond_var` and `cond_shuff` options.
%%%
%%%   3. **Core Function Evaluation**: The `corefunc` is repeatedly called with the shuffled
%%%      inputs to compute the desired output, generating the null distribution. If conditional
%%%      shuffling is enabled, two distributions are generated: one conditioned on the variables
%%%      and one unconditioned.
%%%
%%%   4. **Parallel Processing**: If the `parallel` option is enabled, the function uses
%%%      parallel computing (e.g., `parfor`) to speed up the generation of null distributions.
%%%
%%% Notes:
%%% This function is typically used for statistical hypothesis testing, allowing users to
%%% generate a distribution of random outputs to compare against the actual observed data.
%%% It can be applied in areas such as neuroscience, where transfer entropy or information
%%% flow between neurons is measured.
%%%
%%% Example:
%%% Suppose we have input time series data from two neurons (A and B) and we want to generate a null
%%% distribution to test the significance of a metric like transfer entropy. We can call the function as:
%%%
%%%   null_distribution = create_NullDistribution({A, B}, {'TE'}, @transferEntropy, opts);
%%%
%%% This will generate the null distribution based on shuffling the input time series and
%%% calculating transfer entropy (or any other desired metric) using the core function.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check Inputs, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
default_opts.bias = 'naive';
default_opts.bin_method = {'none'};
default_opts.n_bins = {3};
default_opts.shuff = 0;
default_opts.n_samples = 100;
default_opts.cond_var = {'C'};
default_opts.cond_shuff = false;
default_opts.shuff_var = {'A','B'};
default_opts.parallel = 0;
default_opts.dim_shuffle = {'Trials'};

if nargin < 1
    msg = 'not enough input arguments. See `help shuff_test` for usage info';
    error('shuff_test:notEnoughInput', msg);
end
opts = varargin{1};
if ~isfield(opts, 'supressWarnings')
    opts.supressWarnings = false;
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


if length(opts.cond_shuff) == 2
    doBoth = true;
    opts.cond_shuff = true;
else
    doBoth = false;
end

if opts.cond_shuff
    nInputs = length(inputs);
    cond_vars = zeros(1,nInputs );
    vars_to_check = arrayfun(@(x) char('A' + (x-1)), 1:nInputs, 'UniformOutput', false);
    cond_vars = ismember(vars_to_check, opts.cond_var);
    cond_vars(end) = ismember('S',  opts.cond_var);
    inputs_b = binning(inputs, opts);
    cond_inputs = inputs_b(cond_vars == 1);
else
    cond_inputs = {};
end

opts_funcCall = opts;
opts_funcCall.compute_nulldist = false;
null_distribution = repmat({0}, opts.n_samples, length(outputs));
null_distribution_uncond = repmat({0}, opts.n_samples, length(outputs));
if ~opts.parallel
    for shuffIdx = 1:opts.n_samples
        inputs_sh = shuffle(inputs, opts,cond_inputs);
        values = feval(corefunc, inputs_sh, outputs, opts_funcCall);
        for outIdx = 1:length(outputs)
            null_distribution{shuffIdx, outIdx} = values{outIdx};
        end
    end
    if doBoth
        opts.cond_shuff = false;
        for shuffIdx = 1:opts.n_samples
            inputs_sh = shuffle(inputs, opts,cond_inputs);
            values = feval(corefunc, inputs_sh, outputs, opts_funcCall);
            for outIdx = 1:length(outputs)
                null_distribution_uncond{shuffIdx, outIdx} = values{outIdx};
            end
        end
        null_distribution = {null_distribution, null_distribution_uncond};
    end 
else
    % parfor shuffIdx = 1:opts.n_samples
    %     inputs_sh = shuffle(inputs, opts,cond_inputs);
    %     values = feval(corefunc, inputs_sh, outputs, opts_funcCall);
    %     for outIdx = 1:length(outputs)
    %         null_distribution{shuffIdx, outIdx} = values{outIdx};
    %     end
    % end
    % if doBoth
    %     opts.cond_shuff = false;
    %     parfor shuffIdx = 1:opts.n_samples
    %         inputs_sh = shuffle(inputs, opts,cond_inputs);
    %         values = feval(corefunc, inputs_sh, outputs, opts_funcCall);
    %         for outIdx = 1:length(outputs)
    %             null_distribution_uncond{shuffIdx, outIdx} = values{outIdx};
    %         end
    %     end
    %     null_distribution = {null_distribution, null_distribution_uncond};
    % end

end
end

