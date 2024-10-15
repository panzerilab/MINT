function [outputs, opts] = check_inputs(corefunc,inputs,varargin)


defaultOpts.computeNulldist = false;
defaultOpts.bias = 'naive';
defaultOpts.bin_method = {'none'};
defaultOpts.n_bins = {3};

defaultOutputs_MI      = {'I(A;B)'};
defaultOutputs_cMI     = {'I(A;B|C)'};
defaultOutputs_TE      = {'TE(A->B)'};
defaultOutputs_cTE     = {'TE(A->B|C)'};
defaultOutputs_II      = {'II(A,B,C)'};
defaultOutputs_PID     = {'PID_atoms'};
defaultOutputs_FIT     = {'FIT(A->B;C)'};
defaultOutputs_cFIT    = {'FIT(A->B;C|D)'};
defaultOutputs_H       = {'H(A)', 'H(A|B)'};

if strcmp(corefunc, 'MI')|| strcmp(corefunc, 'cMI') || strcmp(corefunc, 'H')

elseif strcmp(corefunc, 'TE') || strcmp(corefunc, 'cTE')
    defaultOpts.tau = {1};
    defaultOpts.tpres = size(inputs{1}, 2);
    defaultOpts.singleTimepoint = false;
elseif strcmp(corefunc, 'FIT') || strcmp(corefunc, 'cFIT')
    defaultOpts.redundancy_measure = 'I_min';
    defaultOpts.pid_constrained = false;
    defaultOpts.tau = {1};
    defaultOpts.tpres = size(inputs{1}, 2);
    defaultOpts.recall = false; 
elseif strcmp(corefunc, 'II')
    defaultOpts.function = @pidBROJA;
elseif strcmp(corefunc, 'PID')
    defaultOpts.redundancy_measure = 'I_BROJA';
    defaultOpts.pid_constrained = true;
elseif strcmp(corefunc, 'correction')
    defaultOpts.shuff = 30;
    defaultOpts.xtrp = 10;
    defaultOpts.pid_constrained = false;
 
end

if isempty(varargin)
    opts = defaultOpts;
    opts.isBinned = false;
    opts.isChecked = false;
    opts.parallel = false;
    opts.supressWarnings = false;
    outputs = eval(['defaultOutputs_' corefunc]);
elseif length(varargin)==1
    if iscell(varargin{1})
        outputs = varargin{1};
        opts = defaultOpts;
        opts.isBinned = false;
        opts.isChecked = false;
        opts.parallel = false;
        opts.supressWarnings = false;
        if ~opts.supressWarnings
            warning("No opts provided in %s. The function will use the default opts.", corefunc)
        end
    elseif  isstruct(varargin{1})
        opts = varargin{1};
        if ~isfield(opts, 'supressWarnings')
            opts.supressWarnings = false;
        end
        if ~isfield(opts, 'isBinned')
            opts.isBinned = false;
        end
        if ~isfield(opts, 'parallel')
            opts.parallel = false;
        end
        if ~isfield(opts, 'isChecked')
            opts.isChecked = false;
        end
        outputs = eval(['defaultOutputs_' corefunc]);
        if ~opts.supressWarnings
            warning("No outputs list provided in %s. The function will compute the default output.", corefunc)
        end
        if ~opts.isChecked
            default_fields= fieldnames(defaultOpts);
            is_field_present = ismember(default_fields, fieldnames(opts));
            missing_fields = default_fields(~is_field_present);
            for i=1:size(missing_fields,1)
                missing_field_name = missing_fields{i};
                opts.(missing_fields{i}) = defaultOpts.(missing_fields{i});
                if ~opts.supressWarnings
                    if iscell(defaultOpts.(missing_field_name))
                        numericValue = cell2mat(defaultOpts.(missing_field_name));
                    elseif isa(defaultOpts.(missing_field_name), 'function_handle')
                        funcName = func2str(defaultOpts.(missing_field_name));
                        numericValue = funcName;
                    else
                        numericValue = defaultOpts.(missing_field_name);
                    end
                    fprintf('Option "%s" was not specified. Using default value %s.\n', missing_field_name, mat2str(numericValue));
                end
            end
        end
    end
else
    outputs =  varargin{1};
    opts = varargin{2};
    if ~isfield(opts, 'supressWarnings')
        opts.supressWarnings = false;       
    end
    if ~isfield(opts, 'isBinned')
        opts.isBinned = false;             
    end
    if ~isfield(opts, 'parallel')
            opts.parallel = false;
    end
    if ~isfield(opts, 'isChecked')
        opts.isChecked = false;              
    end
    if ~opts.isChecked
        default_fields= fieldnames(defaultOpts);
        is_field_present = ismember(default_fields, fieldnames(opts));
        missing_fields = default_fields(~is_field_present);
        for i=1:size(missing_fields,1)
            missing_field_name = missing_fields{i};
            opts.(missing_fields{i}) = defaultOpts.(missing_fields{i});
            if ~opts.supressWarnings
                if iscell(defaultOpts.(missing_field_name))
                    numericValue = cell2mat(defaultOpts.(missing_field_name));
                elseif isa(defaultOpts.(missing_field_name), 'function_handle')
                    funcName = func2str(defaultOpts.(missing_field_name));
                    numericValue = funcName;
                else
                    numericValue = defaultOpts.(missing_field_name);
                end
                fprintf('Option "%s" was not specified. Using default value %s.\n', missing_field_name, mat2str(numericValue));
            end
        end
    end
end




lastDimLength = [];
nVars = length(inputs);
for var = 1:nVars
    input_var = inputs{var};
    if any(isnan(input_var))
       error('NaN value found in input variable %d in function %s', var, corefunc);
    end
    currentDimLength = size(input_var, ndims(input_var));
    lastDimLength(var) = currentDimLength;
end 
if ~all(lastDimLength == lastDimLength(1))
    error('Input variables have different number of trials in function %s', corefunc);
end
opts.isChecked = true;
end

