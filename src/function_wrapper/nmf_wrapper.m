function output_values = nmf_wrapper(input, varargin)
% *function output_values = nnmf_wrapper(input, outputs, opts)*
%
% Inputs:
%   - input: A cell array containing:
%       - input{1}: Data matrix `a`, with dimensions N-by-M.
%   - varargin: Optional arguments, including:
%       - outputs: A cell array specifying which outputs to return (default is {'all'}).
%       - opts: A structure containing options, including:
%           - newDim: Desired rank for factorization (The target number of reduced dimensions)(k)
%           - algorithm: 'als' (default) for alternating least squares or 'mult' for multiplicative update.
%           - w0: Initial N-by-K matrix for W.
%           - h0: Initial K-by-M matrix for H.
%           - replicates: Number of times to repeat the factorization with different starting points.
%           - options: STATSET structure for additional configuration, e.g., MaxIter, Display, TolFun, TolX.
%
% Outputs:
%   - output_values: A cell array containing the requested output values, which may include:
%       - 'wbest': Optimal factor matrix W (N-by-K).
%       - 'hbest': Optimal factor matrix H (K-by-M).
%       - 'normbest': Final residual norm of the factorization.
%
% Example usage:
%   output_values = nnmf_wrapper({A}, 3, {'wbest', 'hbest'}, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check input, OutputList, and set default options                      %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define default options
defaultOpts.algorithm = 'als';
defaultOpts.replicates = 1;
defaultOpts.w0 = [];
defaultOpts.h0 = [];
defaultOpts.NaN_handling = 'error';

defaultOptionOpts.MaxIter = 100;
defaultOptionOpts.Display = 'off';
defaultOptionOpts.TolFun = 1e-4;
defaultOptionOpts.TolX = 1e-4;
defaultOptionOpts.UseParallel = false;
defaultOptionOpts.UseSubstreams = false;
defaultOptionOpts.Streams = [];

% Check number of arguments and assign variables accordingly
if nargin < 1
    msg = 'Not enough input arguments. See `help nmf_wrapper` for usage info';
    error('NMF:notEnoughInput', msg);
elseif nargin == 1
    opts = defaultOpts;
    outputs = {'all'};
elseif nargin == 2
    if iscell(varargin{1})
        outputs = varargin{1};
        opts = defaultOpts;
    elseif isstruct(varargin{1})
        opts = varargin{1};
        outputs = {'all'};
    end
elseif nargin == 3
    opts = varargin{2};
    outputs = varargin{1};
end
fields = fieldnames(defaultOpts);
for i = 1:numel(fields)
    if ~isfield(opts, fields{i})
        opts.(fields{i}) = defaultOpts.(fields{i});
    end
end
if ~isfield(opts, 'NaN_handling')
    opts.NaN_handling = 'error';
end
input = nan_method(input, opts.NaN_handling);
data = input{1};

if ~isfield(opts,'options')
    opts.options = defaultOptionOpts;
else
    default_fields= fieldnames(defaultOptionOpts);
    is_field_present = ismember(default_fields, fieldnames(opts.options));
    missing_fields = default_fields(~is_field_present);
    for i=1:size(missing_fields,1)
        missing_field_name = missing_fields{i};
        opts.options.(missing_fields{i}) = defaultOptionOpts.(missing_fields{i});
    end
end


if ~isfield(opts, 'newDim')
   error('NNMF:notEnoughInput', ['Not enough input arguments. You need to specify the desired rank for factorization ' ...
       '(The target number of reduced dimensions)(k). See help nmf_wrapper for usage info.']); 
end 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Step 2: Perform NNMF                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Call nnmf with specified options
[wbest, hbest, normbest] = nnmf(data, opts.newDim, ...
    'algorithm', opts.algorithm, ...
    'w0', opts.w0, ...
    'h0', opts.h0, ...
    'replicates', opts.replicates, ...
    'options', opts.options);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Step 4: Outputs                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Define possible outputs and process the output list
possibleOutputs = {'wbest', 'hbest', 'normbest'};
if ismember('all', outputs)
    outputs = possibleOutputs;
end
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    error('NNMF:invalidOutput', 'Invalid Outputs: %s', strjoin(nonMembers, ', '));
end

% Prepare the output
output_values = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'wbest'
            output_values{i} = wbest;
        case 'hbest'
            output_values{i} = hbest;
        case 'normbest'
            output_values{i} = normbest;
    end
end

end




