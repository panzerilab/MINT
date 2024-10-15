function [PID_values, PID_naive, PID_shuff_all] = PID(inputs, varargin)
%%% *function [PID_v_uncorrected, PID_v, p_distr, q_distr, pind_distr, shuff_v] = PID(X, Y,varargin)*
%%%
%%% ### Description
%%% This function computes the partial information decomposition of X1 and
%%% X2 about Y
%%%
%%% ### Inputs:
%%% - *X*: must be a cell array of *nSources X nDimensions X nTrials* elements representing the sources signal.
%%%        nSources is the number of sources (at least two), nDimensions is the number of dimensions of each source,
%%%        and nTrials is the number of trials. The number of dimensions can be different for each source,
%%%        but the number of trials must be the same between sources and also respect to the target.
%%% - *Y*: must be an array of *nDimensions X nTrials* elements representing the target signal
%%% 
%%% - *opts*: options used to calculate PID (see further notes).
%%%
%%% ### Outputs:
%%% - *PID_v*: double array of the corrected information atoms, in the order
%%%   [Shared Information, Unique Information X1, Unique Information X2, Complementary Information].
%%% - *PID_v_uncorrected*: double array of the uncorrected information atoms, in the same
%%%   order as PID_v.
%%% - *p_distr*: The probability distribution of the joint system formed by Y and X.
%%% - *q_distr*: The probability distributionused to calculate the PID terms.
%%% - *pind_distr*: The probability distribution assuming conditional independence of the sources X1 and X2 in respect to Y.
%%% - *shuff_v*: array of size opts.shuff X 4 elements, each row corresponding
%%%   to a single bootstrap iteration and the columns in the same order as PID_v.
%%%
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | Field name and default value        | Description                                                                                           | Possible values
%%% |-------------------------------------|-----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
%%% | opts.bias = `'naive'`               | bias correction method                                                                                | 'naive', 'qe' or 'le'                                                                                                                                        |
%%% | opts.bin_methodX = `'eqpop'`        | binning method for the X signal                                                                       | 'none', 'eqpop', 'eqspace', 'ceqspace', 'gseqspace'                                                                                                          |
%%% | opts.bin_methodY = `'none'`         | binning method for the Y signal                                                                       | 'none', 'eqpop', 'eqspace', 'ceqspace', 'gseqspace'                                                                                                          |
%%% | opts.bin_nbinsX = 3                 | number of bins to reduce the  dimensionality of X                                                          | int > 1                                                                                                                                                      |
%%% | opts.bin_nbinsY = N/A               | number of bins to reduce the dimensionality of Y                                                      | int > 1                                                                                                                                                      |
%%% | opts.shuff = 0                       | number of shuffle operations to perform for significance testing                                    | int >= 0                                                                                                                                                     |
%%% | opts.shuff_variables = N/A           | list of variables to be bootstrapped, specified as a cell array of strings.                           | cell array containing one (or more) strings (`"X1"`, `"X2"` or `"Y"`)                                                                                        |
%%% | opts.shuff_type = N/A                | type of bootstrapping to be applied for each variable                                                 | cell array of same size of opts.shuff_variables, each element contains the shuffling type (either 'all', 'X1conditioned', 'X2conditioned', or 'Yconditioned') |
%%% | opts.old_output = 0                 | flag to indicate if the output will follow the old output from the NIT toolbox (1 if yes)             | int 0 or 1                                                                                                                                                   |
%%% | opts.xtrp = 5                       | how many iterations repetitions of the extrapolation procedure should be performed.                   | int > 0                                                                                                                                                      |
%%% | opts.redundancy_measure = 'I_BROJA' | name of the measure of the redundancy between sources                                                 | 'I_BROJA', 'I_min' or 'I_MMI'                                                                                                                                |
%%% | opts.function = @pidBROJA           | function handle of the redundancy measure between sources, not necessary if redundancy_measure is set | any function handle that accepts a probability distribution as input                                                                                         |
%%% | opts.parallel = 0                   | indicates if the extrapolation and bootstrap procedures run in parallel (1 if yes)                    | int 0 or 1                                                                                                                                                   |
%%%
%%% Bias
%%% - 'naive' (no bias correction)
%%% - 'le' (linear extrapolation)
%%% - 'qe'(quadratic exrapolation)
%%%
%%% Binning methods
%%% - 'none' (no binning)
%%% - 'eqpop' (evenly populated binning)
%%% - 'eqspace' (evenly spaced binning)
%%% - 'ceqspace' (centered evenly spaced binning)
%%% - 'gseqspace' (gaussian evenly spaced binning)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check Inputs, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if length(varargin) > 1
    opts = varargin{2};
    if isfield(opts, 'isChecked')
        if opts.isChecked
            outputs = varargin{1};
        end
    else
        [outputs, opts] = check_inputs('PID',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('PID',inputs,varargin{:});
end

possibleOutputs = {'', 'Syn',  'Red', 'Unq1', 'Unq2', ...
    'Unq', 'PID_atoms', 'Joint', 'Union', 'q_dist'};
if ismember('', outputs)
    outputs = {'Syn', 'Red', 'Unq1', 'Unq2'};
elseif ismember('PID_atoms', outputs)
    outputs(ismember(outputs, 'PID_atoms')) = [];  
    outputs = [outputs, {'Syn', 'Red', 'Unq1', 'Unq2'}];
elseif ismember('all', outputs)
    outputs = {'Syn', 'Red', 'Unq1', 'Unq2', 'Unq', 'Joint', 'Union', 'q_dist'};
end
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('PID:invalidOutput', msg);
end

nSources = length(inputs)-1;
if nSources < 2
    msg = 'At least two sources are required.';
    error('PID:NotEnoughSources', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Binning, reduce dimensions if necessary                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~opts.isBinned == false
    inputs_b = binning(inputs, opts);   
    opts.isBinned = true;
else
    inputs_b = inputs;
end

inputs_1d = inputs_b;
nTrials_comp = size(inputs_1d{1},2);
for var = 1:nSources+1
   sizeVar = size(inputs_1d{var});
   nTrials = sizeVar(end);
   if nTrials ~= nTrials_comp
       msg = 'Inconsistent input size. Number of Trials must be consistent.';
       error('PID:Invalid Input', msg);
   end 
   if sizeVar(1) > 1
       inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
   end 
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 3.A: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
corr = opts.bias;
corefunc = @PID;
if ~strcmp(corr, 'naive')
    [PID_values, PID_naive, PID_shuff_all] = correction(inputs_1d, outputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 3.B: Compute required PID Atoms                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

switch opts.redundancy_measure
    case 'I_BROJA'
        if nSources == 2
            opts.function = @pidBROJA;
        else
            warning('BROJA can be used for two sources only, switching to I_min')
            opts.redundancy_measure = 'I_min';
            opts.function = @pidimin;
        end
    case 'I_MMI'
        opts.function = @pidimmi;
    case 'I_min'
        opts.function = @pidimin;
end
opts.multidim = true;
sourceInput = cat(1, inputs_1d{1:end-1});
[p_distr] = prob_estimator({inputs_1d{end}, sourceInput}, {'P(A,B,C)'}, opts);
if strcmp(opts.redundancy_measure,'I_BROJA')
    [PID_terms, table_prob] = opts.function(p_distr);
else
    PID_terms = opts.function(p_distr);
end

if opts.pid_constrained
    I1 = cell2mat(MI({inputs_1d{1}, inputs_1d{end}}, {'I(A;B)'}, opts));
    I2  = cell2mat(MI({inputs_1d{2}, inputs_1d{end}}, {'I(A;B)'}, opts));
    I12 = cell2mat(MI({cat(1, inputs_1d{1}, inputs_1d{2}), inputs_1d{end}}, {'I(A;B)'}, opts));
    if ~isfield(opts, 'chosen_atom')
        opts.chosen_atom = 'red';
    end
    if strcmp(opts.chosen_atom, 'red')
        red = PID_terms(1);
        PID_terms = [red, I1-red, I2-red, I12-I1-I2+red];
    elseif strcmp(opts.chosen_atom, 'un1')
        un1 = PID_terms(2);
        PID_terms = [I1-un1, un1, I2-I1+un1, I12-I2-un1];
    elseif strcmp(opts.chosen_atom, 'un2')
        un2 = PID_terms(3);
        PID_terms = [I2-un2, I1-I2+un2, un2, I12-I1-un2];
    elseif strcmp(opts.chosen_atom, 'syn')
        syn = PID_terms(4);
        PID_terms = [I1+I2-I12+syn, I12-I2-syn, I12-I1-syn , syn];
    end
end


for i = 1:length(outputs) 
    output = outputs{i};
    switch output
        case 'Syn'
            PID_values{i} = PID_terms(4);
        case 'Red'
            PID_values{i} = PID_terms(1);
        case 'Unq1'
           PID_values{i} = PID_terms(2);
        case 'Unq2'
            PID_values{i} = PID_terms(3);
        case 'Unq'
            PID_values{i} = PID_terms(2)+PID_terms(3);
        case 'Joint'
            PID_values{i} = sum(PID_terms);
        case 'Union'
            PID_values{i} = sum(PID_terms)-PID_terms(4);
        case 'q_dist'
            if strcmp(opts.bias,'naive') && strcmp(opts.redundancy_measure,'I_BROJA') && istable(table_prob)
                q_distr = zeros(size(p_distr));
                for irow=1:height(table_prob)
                    x = table_prob{irow,'X'};
                    y = table_prob{irow,'Y'};
                    z = table_prob{irow,'Z'};
                    q_distr(x, y, z) = table_prob{irow, 'qpos'};
                end
                PID_values{i} = q_distr;
            else 
                PID_values{i} = NaN;
                if ~opts.supressWarnings
                    warning("q_dist can only be computed when the redundancy measure is set to 'I_Broja'.")
                end 
            end 
    end
end
PID_naive = PID_values;
end 
