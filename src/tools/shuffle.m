function  inputs_sh = shuffle(inputs, varargin)
% Shuffle function for neural data with optional consistency constraints.
% This function shuffles neural data based on specified
% conditions, allowing customization through input parameters.

% Input:
%   - cond_input: Behavioral data used for consistency constraints
%   - neural_data: nTrials x nObjects x ..., neural data to be shuffled
%   - consistency: Binary flag (0 or 1) indicating whether to apply
%                  consistency constraints for shuffling across trials (0 for no constraints)
%   - index: Binary vector indicating which dimensions to shuffle

% Output:
%   - shuffled_data: Shuffled neural data based on specified conditions

% Description:
% The function shuffles neural data based on the provided input parameters.
% - If consistency is set to 1, the shuffling is consistent across trials,
%   ensuring that unique values in the behavioral data remain grouped together.
% - The index vector indicates which dimensions should be shuffled.
%   index(2) == 1 signifies shuffling over the second dimension (neuralObjects),
%   potentially altering the probability distributions of individual objects.

default_opts.cond_shuff = false;
default_opts.shuff_var = {'A','B'}; % A B C D S All
default_opts.dim_shuffle = {'Trials'}; %Trials Timepoints, Objects or number
if nargin < 1
    msg = 'not enough input arguments. See `help shuff_test` for usage info';
    error('shuffle:notEnoughInput', msg);
end
opts = varargin{1};
if nargin > 2
    cond_input = varargin{2};
end
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
% Create Array for shuffle Vars and Cond Vars
nVars = length(inputs);

vars_to_check = arrayfun(@(x) char('A' + (x-1)), 1:nVars, 'UniformOutput', false);
shuffle_vars = ismember(vars_to_check, opts.shuff_var);
shuffle_vars(end) = ismember('S', opts.shuff_var);

% Check inputs
dimsA = size(inputs{1});
if opts.cond_shuff
    for var = 2:nVars
        dims = size(inputs{var});
        if dims(end) ~= dimsA(end)
            letter_var = char('A' + var - 1);
            error('Number of trials for Var A and Var %s are inconsistent.', letter_var);
        end
    end
end

nDims = length(dimsA);
% Create Array for dimensions to shuffle
shuffle_dim = zeros(1, length(dimsA));
for i = 1:length(opts.dim_shuffle)
    dim_shuffle_entry = opts.dim_shuffle{i};
    if ischar(dim_shuffle_entry)
        if strcmp(dim_shuffle_entry, 'Objects')
            shuffle_dim(1) = 1;
        elseif strcmp(dim_shuffle_entry, 'Timepoints')
            shuffle_dim(2) = 1;
        elseif strcmp(dim_shuffle_entry, 'Trials')
            shuffle_dim(end) = 1;
        end
    end
end
if opts.cond_shuff
    cond_input_data = [];
    for var = 1:length(cond_input)
        cond_input_data =[cond_input_data; cond_input{var}];
    end
end
inputs_sh = inputs;
if nDims > 2
    nTimepoints = dimsA(2);
end

if length(nDims) == 2
    for var = 1:nVars
        if shuffle_vars(var) == 1  
            if shuffle_dim(end)==1
                if opts.cond_shuff
                    shuffled_data = (shuffle_core(cond_input',inputs_sh{var}',1, 1)');
                else 
                    shuffled_data = (shuffle_core(0,inputs_sh{var}',0, 1));  
                end 
            end
            inputs_sh{var} = shuffled_data;
        end
    end
else
    for var = 1:nVars
        if shuffle_vars(var) == 1
            inputs_data = inputs_sh{var};
            dims = size(inputs_data);
            shuffled_data = inputs_sh{var};
            for tP = 1:nTimepoints
                data_tP = squeeze(inputs_data(:,tP,:));
                data_tP = reshape(data_tP, dims(3), dims(1));
                if opts.cond_shuff
                    if length(size(cond_input_data))>2
                        cond_input_tp = squeeze(cond_input_data(:,tP,:));
                        cond_input_tp = reshape(cond_input_tp, size(cond_input,3), size(cond_input,1));
                    else 
                        cond_input_tp = cond_input_data;
                    end
                end 
                if shuffle_dim(end)==1
                    if opts.cond_shuff                       
                        shuffled_data_tmp = (shuffle_core(cond_input_tp', data_tP,1, 1))';
                    else                           
                        shuffled_data_tmp = (shuffle_core(0, data_tP,0, 1))';  
                    end 
                end
                shuffled_data(:,tP,:) =  shuffled_data_tmp;
            end 
            inputs_sh{var} = shuffled_data;
       end
    end
end 
end