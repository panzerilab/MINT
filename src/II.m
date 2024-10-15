function [II_values, II_naive, II_nullDist, atom1, atom2] = II(inputs, varargin)
%%% *function [II, PID_SR_C, PID_RC_S] = II(S, R, C, varargin)*
%%%
%%% ### Description
%%% This function computes the intersection information (II), either
%%% naive or bias-corrected, based on the input data.
%%%
%%% ### Inputs:
%%% - *S*: An array of size *nDimsS x nTrials* representing the discrete 
%%%        stimulus presented in each trial.
%%% - *R*: An array of size *nDimsR x nTrials* representing the response
%%%        of each of the *nDimsR* dimensions for each trial.
%%% - *C*: An array of size *nDimsC x nTrials* representing the discrete
%%%        choice made by the subject in each trial.
%%% - *opts*: Optional settings used to calculate II (see further notes).
%%%
%%% ### Outputs:
%%% - *II_v*: final value of intersection information. If opts.bias is 'naive', it's the
%%% uncorrected value. If opts.bias is anything else, it is the value after
%%% bias correction.
%%% - *II_naive*: uncorrected value of intersection information. 
%%% - *shuff_v*: intersection information values for each permutation.
%%%
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% - *opts.bias*: Specifies the bias correction method.
%%%   Allowed values:
%%%     - `'naive'`: No bias correction (default)
%%%     - `'le'`: Linear extrapolation
%%%     - `'qe'`: Quadratic extrapolation
%%%
%%% - *opts.bin_methodS*, *opts.bin_methodR*, *opts.bin_methodC*:
%%%   Specify the binning method for the stimulus (*S*), response (*R*), 
%%%   and choice (*C*), respectively. The following binning methods 
%%%   are allowed:
%%%     - `'none'`: No binning (default for *S* and *C*)
%%%     - `'eqpop'`: Evenly populated binning (default for *R*)
%%%     - `'eqspace'`: Evenly spaced binning
%%%   See the documentation for the `binr` function for more details.
%%%
%%% - *opts.n_binsS*: Specifies the number of bins used to reduce the 
%%%   dimensionality of the stimulus. Must be an integer > 1 (default: 3).
%%%
%%% - *opts.n_binsR*: Specifies the number of bins used to reduce the 
%%%   dimensionality of the response. Must be an integer > 1 (default: 3).
%%%
%%% - *opts.n_binsC*: Specifies the number of bins used to reduce the 
%%%   dimensionality of the choice. Must be an integer > 1 (default: 3).
%%%
%%% - *opts.shuff*: Specifies the number of bootstrap operations for 
%%%   significance testing. The bootstraps are performed independently on
%%%   each variable listed in `opts.shuff_variables`. Must be an integer ≥ 0
%%%   (default: 0).
%%%
%%% - *opts.shuff_variables*: A cell array of strings specifying which
%%%   variables to bootstrap. Can include one or more of the following:
%%%     - `"S"`: Stimulus
%%%     - `"R"`: Response
%%%     - `"C"`: Choice
%%%   If multiple bootstraps are requested, all specified variables are 
%%%   shuffled.
%%%
%%% - *opts.shuff_type*: Specifies the type of bootstrapping to apply to 
%%%   each variable in `opts.shuff_variables`. Each element of the cell array 
%%%   should specify the type of shuffling:
%%%     - `'all'`: Shuffles all values of the variable across trials.
%%%     - `$VAR$conditioned`: Shuffles the trials conditioned on the values 
%%%       of the variable `$VAR$`, where `$VAR$` can be `"S"`, `"R"`, or `"C"`.
%%%
%%% - *opts.shuff_bias*: Specifies the bias correction method to be applied 
%%%   for each variable in `opts.shuff_variables`. Each element should 
%%%   specify one of the bias correction methods described in `opts.bias`.


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
        [outputs, opts] = check_inputs('II',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('II',inputs,varargin{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%               Step 2: Binning, reduce dimensions if necessary                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if ~opts.isBinned
    inputs_b = binning(inputs ,opts);
    opts.isBinned = true;
else
    inputs_b = inputs;
end
inputs_1d = inputs_b;

for var = 1:length(inputs_b)
    sizeVar = size(inputs_1d{var});
    if sizeVar(1) > 1
        inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
    end
end


possibleOutputs = {'II(B,A,C)','II(A,B,C)', 'II(A,C,B)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('II:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%            Step 3.A: Bias correction and nullDist if requested                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
II_nullDist = 0;
corr = opts.bias;
corefunc = @II;
if any(opts.computeNulldist)
        II_nullDist = create_NullDistribution(inputs, outputs, @II, nullDist_opts);
end
if ~strcmp(corr, 'naive')  
    opts.compute_nulldist = false;
    [II_values, II_naive] = correction(inputs_1d, outputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Step 3.B: Compute Probability Distributions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_distributions = struct( ...
    'II_BAC', {{'P(B,A,C)'}}, ...
    'II_ABC', {{'P(A,B,C)'}}, ...
    'II_ACB', {{'P(A,C,B)'}} ...
    );

opts.multidim = true;

required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'II(B,A,C)'
            required_distributions = [required_distributions, entropy_distributions.II_BAC{:}];
        case 'II(A,B,C)'
            required_distributions = [required_distributions, entropy_distributions.II_ABC{:}];
        case 'II(A,C,B)'
            required_distributions = [required_distributions, entropy_distributions.II_ACB{:}];
    end
end

prob_dists = {};
for i = 1:length(required_distributions)
    switch required_distributions{i}
        case 'P(B,A,C)'
            prob_dists{i} =  prob_estimator({inputs_1d{2},inputs_1d{1},inputs_1d{3}}, {'P(A,B,C)'}, opts);
        case 'P(A,B,C)'
            prob_dists{i} =  prob_estimator({inputs_1d{1},inputs_1d{2},inputs_1d{3}}, {'P(A,B,C)'}, opts);
        case 'P(A,C,B)'
            prob_dists{i} =  prob_estimator({inputs_1d{1},inputs_1d{3},inputs_1d{2}}, {'P(A,B,C)'}, opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 3.C: Compute requested Output Values                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the TE values
II_values = cell(1, length(outputs));
atom1 = cell(1, length(outputs));
atom2 = cell(1, length(outputs));

for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'II(B,A,C)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(B,A,C)')};
            atoms = ii_core(Prob_d, opts);
            II_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
        case 'II(A,B,C)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(A,B,C)')};
            atoms = ii_core(Prob_d, opts);
            II_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
        case 'II(A,C,B)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(A,C,B)')};
            atoms = ii_core(Prob_d, opts);
            II_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
    end
end
II_naive = II_values;
end