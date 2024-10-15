function [inputs_b, edg] = binning(inputs, varargin)
%%% *function [inputs_b] = binning(inputs, varargin)*
%%%
%%% The binning function takes a set of input data arrays and applies various binning methods 
%%% to discretize continuous values into specified categories or bins. This is useful for 
%%% data analysis, allowing continuous data to be represented in a more manageable form.
%%%
%%% Inputs:
%%%   - inputs: A cell array containing the input data arrays to be binned. Each cell can 
%%%             represent a different variable or data structure. The function can handle 
%%%             multiple variables simultaneously.
%%%
%%%   - varargin: Optional arguments, passed as a structure. Fields may include:
%%%              - bin_method: A cell array specifying the binning method for each input. 
%%%                            Possible methods include:
%%%                            'none'       : No binning applied (default).
%%%                            'eqspace'    : Equal space binning.
%%%                            'threshold'  : Binning based on specified thresholds.
%%%                            'eqpop'      : Equal population binning.
%%%              
%%%              - n_bins: A cell array specifying the number of bins to use for binning 
%%%                        for each input. It can be a single integer or a cell array.
%%%                        Default is {3}.
%%%
%%%              - computeOptimBins: Boolean (true/false) indicating whether to compute 
%%%                                  optimal bins automatically (default is false).
%%%              
%%%              - computeBinsMethod: Specifies the method to compute optimal bins, 
%%%                                   either 'Scott' or 'Freedman' (default is 'Scott' 
%%%                                   if computeOptimBins is true).
%%%
%%%              - supressWarnings: Boolean (true/false) to suppress warning messages. 
%%%                                 Default is false, meaning warnings will be shown.
%%%
%%% Outputs:
%%%   - inputs_b: A cell array containing the binned data arrays corresponding to the 
%%%                specified binning methods and parameters.
%%%
%%% Note: 
%%% The binning method can be tailored for each variable, allowing for flexible data processing 
%%% strategies.
%%%
%%% Custom Binning Strategies:
%%% Users can define their own custom binning strategies by creating a separate MATLAB function 
%%% that implements the desired binning logic. This function should take the following inputs:
%%%   - inputs: The data to be binned.
%%%   - nbins: The number of bins to apply.
%%%   - opts: Options and parameters for binning.
%%% The custom function should output the binned data. The custom binning strategy can then 
%%% be invoked by specifying the function name in the `bin_method` field. The function must 
%%% be placed in the `tools` directory and follow the naming convention `your_custom_method.m`.
%%%
%%% EXAMPLE 1
%%% Suppose we have a set of continuous data arrays representing different variables:
%%% inputs = {randn(1, 100), rand(1, 100)}; 
%%% To discretize the data into equal population bins with 5 bins for each variable, 
%%% the function can be called as:
%%%
%%% opts.bin_method = {'eqpop', 'eqpop'};
%%% opts.n_bins = {5, 5};
%%% inputs_b = binning(inputs, opts);
%%%
%%% EXAMPLE 2
%%% inputs = {randn(2, 100), rand(1, 100)};
%%% To discretize the data into equal population bins with:
%%%  - 5 bins for the first variable, first dimension 
%%%  - 3 bins for the first variable, second dimension 
%%%  - 5 bins for the second variable
%%% the function can be called as:
%%% opts.bin_method = {'eqpop', 'eqpop'};
%%% opts.n_bins = {[5,3], 5};
%%% inputs_b = binning(inputs, opts);
%%% The binned data for each variable will be returned in the inputs_b cell array.


if nargin < 1
    error("binning.m: not enough input arguments.");
end

default_opts.bin_method = 'none';
default_opts.n_bins = {3};
default_opts.supressWarnings = false;


if nargin == 1
    opts = default_opts;
    msg = "No opts specified, there will be no binning except unique values is bigger than 100.";
    warning('binning.m:undefined opts.',msg);
else
    opts = varargin{1};
    if ~isfield(opts, 'supressWarnings')
        opts.supressWarnings = false;
    end
    default_fields= fieldnames(default_opts);
    is_field_present = ismember(default_fields, fieldnames(opts));
    missing_fields = default_fields(~is_field_present);
    for i=1:length(missing_fields)
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

end

if ~isfield(opts, 'computeOptimBins')
    opts.computeOptimBins = false;
end
if opts.computeOptimBins == true && ~isfield(opts, 'computeBinsMethod')
    opts.computeBinsMethod = 'Scott';
end

nVars = length(inputs);
if length(opts.bin_method) < nVars
    opts.bin_method((end+1):nVars) = repmat(opts.bin_method(end), 1, nVars - length(opts.bin_method));
end
if length(opts.n_bins) < nVars && opts.computeOptimBins == false
     opts.n_bins((end+1):nVars) = repmat(opts.n_bins(end), 1, nVars - length(opts.n_bins));
end
inputs_b = inputs;


for var = 1:nVars

    unique_vals = unique(inputs{var})';

    if strcmp(opts.bin_method{var}, 'none') && size(unique_vals,2) > 100
        opts.bin_method{var} = 'eqpop';
        opts.n_bins{var} = 100;
    end

    if ~strcmp(opts.bin_method{var}, 'none')
        data = inputs{var};
        [nDims, ~] = size(data);

        if opts.computeOptimBins == true
            switch opts.computeBinsMethod
                case 'Scott'
                    nBinsFunc = @scottsRule;
                case 'Freedman'
                    nBinsFunc = @freedmanDiaconisRule;
                otherwise
                    error('Invalid method to compute the optimal number of bins specified. Options are ´Scott´ or ´Freedman´');
            end

            for dim = 1:nDims
                if length(size(data))== 2
                    nbins_tmp(dim) = nBinsFunc(data(dim,:));
                else
                     nbins_tmp(dim) = nBinsFunc(data(dim,:,:));
                end
            end
            opts.n_bins{var} = nbins_tmp;
        end

        nbins = opts.n_bins{var};
        bin_method = opts.bin_method{var};


        nbins_length = length(nbins);

        if nbins_length ~= nDims
            if nbins_length == 1
                nbins = repmat(nbins, nDims, 1);  % Verwende repmat für eine bessere Lesbarkeit
            else
                error('binning.m:N_binsLengthMismatch', ...
                    'The number of bins (%d) exceeds the number of dimensions (%d). Please ensure that the number of bins matches the size of the input data.', ...
                    nbins_length, nDims);
            end
        end

        data_binned = inputs{var};
        for dim = 1:nDims
            if length(size(data))== 2
                data_dim = data(dim, :);
            else
                data_dim = data(dim, :,:);
            end
            switch bin_method
                case 'eqspace'
                    minValue = min(data_dim(:));
                    maxValue = max(data_dim(:));

                    if ~isfield(opts, 'range')
                        leftEdg =  minValue;
                        rightEdg =  maxValue;
                    else
                        range = opts.range;
                        if length(range)==2 && range(1)<range(2)
                            lftEdg = range(1);
                            rgtEdg = range(2);
                        else
                            error('binning.m_eqspace:InvalidRange. Invalid range specified for var (%d) in dimension (%d).', var, dim);
                        end

                        if minValue<leftEdg || maxValue>rightEdg
                            error('binning.m_eqspace:InvalidRange. Values out of provided range for var (%d) in dimension (%d).', var, dim);
                        end
                    end
                    edg = linspace(leftEdg, rightEdg, nbins(dim)+1);
                    edg(nbins(dim)+1) = edg(nbins(dim)+1) + 1;
                case 'threshold'
                    if ~isfield(opts, 'threshold')
                        error('binning.m_threshold:MissingTreshold. Thresholds must be specified for var (%d) in dimension (%d).', var, dim);
                    end
                    tresholds = opts.threshold;
                    if ~isvector(tresholds) || length(tresholds) < 2
                        error('binning.m_treshold:InvalidTresholds. At least two tresholds must be specified for var (%d) in dimension (%d).', var, dim);
                    end
                    minValue = min(data_dim(:));
                    maxValue = max(data_dim(:));
                    if any(tresholds < minValue) || any(tresholds > maxValue)
                        if ~opts.supressWarnings
                            warning('binning.m_treshold:OutOfRange. Some tresholds are out of the range of data for var (%d) in dimension (%d).', var, dim);
                        end
                    end
                    if length(tresholds) == 1 || isscalar(tresholds) && (isnumeric(tresholds) || isinteger(tresholds))
                        edg = [-inf, tresholds(1), inf];
                    elseif length(tresholds) >= 2
                        edg = [-inf, tresholds(:)', inf];
                    else
                        error('binning.m_treshold:InvalidTresholds. At least one valid threshold must be specified for var (%d) in dimension (%d).', var, dim);
                    end

                case 'eqpop'
                    data_tmp = data(dim, :);
                    unique_vals_dim = unique(data_tmp);
                    N = size(data_tmp,2);
                    binSize = floor(N / nbins(dim));
                    if length(unique_vals_dim)<nbins(dim)
                        error('binning.m_eqpop:Too many bins specified for var (%d) in dimension (%d).', var, dim);
                    end
                    if sum(data_tmp == mode(data_tmp)) > binSize
                        if ~opts.supressWarnings
                            warning('binning.m_eqpop:tooFrequentResponse',"Using equally populated bins with a response array showing at least one response that is more frequent than nSamples/nBins. Therefore the bins are not exactly equipopulated")
                        end
                    end

                    [~ , sortedIndxs] = sort(data_tmp);

                    idxEdges = round(linspace(1, N+1 , nbins(dim) + 1));
                    idxEdges(end) = N;

                    [sortedData, ~] = sort(data_tmp);

                    for i = 2:length(idxEdges)
                        while sortedData(sortedIndxs(idxEdges(i))) == sortedData(sortedIndxs(idxEdges(i - 1)))
                            if idxEdges(i) >= size(data_tmp,2)
                                break;
                            end
                            idxEdges(i:end) = idxEdges(i:end) + 1;
                            idxEdges(idxEdges > size(data_tmp,2)) = size(data_tmp,2);
                        end
                    end

                    idxEdges(idxEdges > size(data_tmp,2)) = size(data_tmp,2);

                    edg = data_tmp(sortedIndxs(idxEdges));

                    edg = unique(edg);

                    idxEdges(idxEdges == length(data_tmp)) = idxEdges(idxEdges == length(data_tmp)) + 1;
                    binCounts = diff(idxEdges);

                    edg = unique(edg);

                    if length(edg) < nbins(dim)+1
                        reducedBins = length(edg)-1;
                        if ~opts.supressWarnings
                            warning(sprintf('binning.m_eqpop: Bin number reduced to %d', reducedBins));
                        end
                    end

                    if any(abs(diff(binCounts)) > 0.05 * max(binCounts))
                        nonZeroStr = '';
                        for i = 1:length(binCounts)
                            nonZeroStr = [nonZeroStr, sprintf('Bin %d: (%d), ', i, binCounts(i))];
                        end

                        nonZeroStr = strtrim(nonZeroStr);
                        if endsWith(nonZeroStr, ',')
                            nonZeroStr = nonZeroStr(1:end-1);
                        end
                        if ~opts.supressWarnings
                            warning('values:unequalSpacing', 'Bin Sizes: %s', nonZeroStr);
                        end
                    end
                otherwise
                    tools_dir = fullfile(pwd, 'tools');
                    func_name = strcat(corefunc, '.m');
                    function_path = fullfile(tools_dir, func_name);
                    if exist(function_path, 'file') == 2
                        % You can change the outputs and inputs of your function here,
                        % but make sure that it fits the output definition of this
                        % function.
                        [data_binned] = feval(func_name, inputs, nbins(dim), opts);
                    else
                        available_functions = {'eqspace','threshold', 'eqpop'};
                        error(['The function "%s" is not defined in the tools folder.\n', ...
                            'Available options are: %s.\n', ...
                            'You can define additional functions by defining your own correction function in a "%s.m" file in the tools folder.'], ...
                            func_name, strjoin(available_functions, ', '), func_name);
                    end
            end
            if length(size(data))== 1
                data_binned(dim) = discretize(data_dim, edg);
            elseif length(size(data))== 2
                data_binned(dim,:) = discretize(data_dim, edg);
            elseif length(size(data))== 3
                data_binned(dim,:,:) = discretize(data_dim, edg);
            end
        end
        inputs_b{var} = data_binned;
    end
end
end
