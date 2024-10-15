function [FIT_values, FIT_naive, FIT_nullDist, atom1, atom2] = FIT(inputs, varargin)
%%% *function [FIT_values, FIT_naive, FIT_shuff_all] = FIT(inputs,outputs, varargin)*
%%%
%%% ### Description
%%% This function computes the Feature-specific (i.e. $S$-specific) Information Transmission (FIT) between a *source* $X$ and a *receiver* $Y$. NIT provides tools for both naive or bias-corrected estimation, given input data.
%%%
%%% ### Inputs:
%%% - *S*: must be an array of *nDimsS X nTrials * elements representing the discrete value of the stimulus presented in each trial.
%%% - *hX*: must be an array of *nDimsX X nTimepoints X nTrials* response matrix describing the past of the emitter variables on each of the *nDims* dimensions for each trial.
%%% - *hY*: must be an array of *nDimsY X nTimepoints X nTrials* response matrix describing the past of the receiver variable on each of the *nDims* dimensions for each trial.
%%% - *Y*: must be an array of *nDimsY X nTimepoints X nTrials* response matrix describing the response of each of the *nDims* dimensions for each trial.
%%% - *opts*: options used to calculate FIT (see further notes).
%%%
%%% ### Outputs:
%%% - *FIT*: data structure containing naive, and bias corrected (if required) estimation of the $S$-related Information Transfer between the *source* $X$ and the *receiver* $Y$. If calculated, the structure contains as well the shuff estimate. In case of multiple bootstrapping options (i.e. in case of `length(opts.shuff_variables) > 1`) multiple bootstrap estimates are returned.
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%

%%% | field                                | description                                                                                                                                                                                                                                             | allowed values                                                                                                                                                                                                                                     | default   |
%%% |--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
%%% | opts.bias                            | specifies the bias correction method                                                                                                                                                                                                                    | `'naive'` (no bias correction)<br>`'le'` (linear extrapolation)<br>`'qe'` (quadratic extrapolation)                                                                                                                                                | `'naive'` |
%%% | opts.bin_methodX                     | specifies the binning method for the X signals                                                                                                                                                                                                          | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_methodY                     | specifies the binning method for the Y signals                                                                                                                                                                                                          | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.bin_methodS                     | specifies the binning method for the stimulus signals                                                                                                                                                                                                   | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details   | `'eqpop'` |
%%% | opts.n_binsX                         | number of bins to be used to reduce the dimensionality of the response X                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                            | 3         |
%%% | opts.n_binsY                         | number of bins to be used to reduce the dimensionality of the response Y                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                            | 3         |
%%% | opts.n_binsS                         | number of bins to be used to reduce the dimensionality of the stimulus S                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                            | 2         |
%%% | opts.shuff                           | number of bootstrap operations to perform for significance testing (those will be performed independently on each of the variables listed in `opts.shuff_variables`)                                                                                                              | int >= 0                                                                                                                                                                                                                                                                   | 0         |
%%% | opts.shuff_variables                 | list of variables to be bootstrapped, specified as a cell array of strings. If multiple bootstrapping operations are requested, `opts.shuff_variables` can contain multiple variables. In this case a `opts.shuff_type` option should be specified for each `opts.shuff_variables`  | cell array containing one (or more) strings (`"S"`, `"R"` or `"C"`) corresponding to all variables to be bootstrapped                                                                                                                                                      | N/A       |
%%% | opts.shuff_type                      | type of bootstrapping to be applied for each variable (`'all'` shuffles all values of the corresponding variable in `opts.shuff_variables` across trials, while `'$VAR$conditioned'` shuffles trials by conditioning on values of the variable specified in the substring `$VAR$`)| cell array of same size of opts.shuff_variables, each element contains the type of shuffling to be performed for the variable (either "all"` or `"Cconditioned"`, `"Rconditioned"`, or `"Sconditioned"`)                                                                    | N/A       |
%%% | opts.shuff_bias                      | bias correction method to be applied for each variable (same values than in opts.bias)
%%% | opts.taux                            | lag considered for the causing signal. Lag is specified as **strictly negative** integer.                                                                                                                                                               | int < 0                                                                                                                                                                                                                                                         | 3         |
%%% | opts.tauy                            | lag considered for the caused signal. Lag is specified as **strictly negative** integer.
%%% | opts.tpres                           | timepoint which should be considered as present in the computation (integer)
%
% | 3         |
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
        [outputs, opts] = check_inputs('FIT',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('FIT',inputs,varargin{:});
end

nullDist_opts = opts;
% nullDist_opts.compute_nulldist = false;

if ~opts.recall
    nSources = length(inputs)-1;
    if nSources < 2
        msg = 'Two sources are required.';
        error('FIT:NotEnoughSources', msg);
    end
    DimsA = size(inputs{1});
    DimsB = size(inputs{2});
    DimsS = size(inputs{3});
    if DimsA(end) ~= DimsB(end) || DimsA(end) ~= DimsS(end)
        msg = sprintf('The number of trials for A (%d), B (%d) and S (%d) are not consistent. Ensure both variables have the same number of trials.',DimsA(3),DimsB(3),DimsS(3));
        error('FIT:InvalidInput', msg);
    end
    nTrials = DimsA(3);

    nTimepointsA = DimsA(2);
    nTimepointsB = DimsB(2);

    if iscell(opts.tau)
        if numel(opts.tau) == 1
            opts.tau{2} = opts.tau{1};
        end
    else
        opts.tau = {opts.tau};
        opts.tau{2} = opts.tau{1};
    end
    if iscell(opts.tpres)
        if numel(opts.tpres) == 1
            opts.tpres{2} = opts.tpres{1};
        end
    else
        opts.tpres = {opts.tpres};
        opts.tpres{2} = opts.tpres{1};
    end

    Atau = opts.tau{1};
    Btau = opts.tau{2};
    Atau = [0; Atau(:)];
    Btau = [0; Btau(:)];


    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %        Step 2: Prepare Data (Shift based on tau,tpres/bin/reduce dims)        %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    A_delayed = zeros(length(Atau), 1, nTrials);
    for tau = 1:length(Atau)
        index = opts.tpres{1} - Atau(tau);
        if index > 0 && index <= nTimepointsA
            A_delayed(tau,:, :) = inputs{1}(:, index, :);
        else
            msg = "The specified delay for A is not available for the timeseries.";
            error('FIT:InvalidInput', msg);
        end
    end
    B_delayed = zeros(length(Btau),1, nTrials);
    for tau = 1:length(Btau)
        index = opts.tpres{2} - Btau(tau);
        if index > 0 && index <= nTimepointsB
            B_delayed(tau,:, :) = inputs{2}(:, index, :);
        else
            msg = "The specified delay for B is not available for the timeseries.";
            error('TE:InvalidInput', msg);
        end
    end

    A_pres = squeeze(A_delayed(1, :, :));
    B_pres = squeeze(B_delayed(1, :, :));
    B_past = squeeze(B_delayed(2:end, :, :));
    A_past = squeeze(A_delayed(2:end, :, :));
    S = inputs{end};
    A_pres = reshape(A_delayed(1, :, :), 1, DimsA(3));  % Shape [1, nTrials]
    B_pres = reshape(B_delayed(1, :, :), 1, DimsB(3));  % Shape [1, nTrials]

    A_past = reshape(A_delayed(2:end, :, :), size(A_delayed, 1)-1, DimsA(3));  % Shape [nDims-1, nTrials]
    B_past = reshape(B_delayed(2:end, :, :), size(B_delayed, 1)-1, DimsB(3));  % Shape [nDims-1, nTrials]

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %               Step 3: Binning, reduce dimensions if necessary                 %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Reshape opts.bin_method
    if numel(opts.bin_method) == 1
        opts.bin_method = repmat(opts.bin_method, 1, 5); % Repeat 5 times
    elseif numel(opts.bin_method) == 2
        opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 3)]; % First 2x, second 3x
    else
        opts.bin_method = [repmat(opts.bin_method(1), 1, 2), repmat(opts.bin_method(2), 1, 2), opts.bin_method(3)];  % 2x first, 2x second, 1x third
    end

    % Reshape opts.n_bins
    if numel(opts.n_bins) == 1
        opts.n_bins = repmat(opts.n_bins, 1, 5); % Repeat 5 times
    elseif numel(opts.n_bins) == 2
        opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 3)]; % First 2x, second 3x
    else
        opts.n_bins = [repmat(opts.n_bins(1), 1, 2), repmat(opts.n_bins(2), 1, 2), opts.n_bins(3)];  % 2x first, 2x second, 1x third
    end

    % Bin the data
    if ~opts.isBinned
        inputs_b = binning({A_pres,A_past,B_pres,B_past, S} ,opts);
        opts.isBinned = true;
    end 
    inputs_1d = inputs_b;
    % Reduce the Dimensions if necessary
    for var = 1:length(inputs_b)
        sizeVar = size(inputs_1d{var});
        if sizeVar(1) > 1
            inputs_1d{var} = reduce_dim(inputs_1d{var}, 1);
        end
    end
     opts.recall = true;
else
    inputs_1d = inputs;
end


possibleOutputs = {'FIT(A->B;S)', 'FIT(B->A;S)'};
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('FIT:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                    Step 4.A: Bias correction if requested                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
FIT_nullDist = 0;
corr = opts.bias;
corefunc = @FIT;
if any(opts.compute_nulldist)
        FIT_nullDist = create_NullDistribution(inputs, outputs, @FIT, nullDist_opts);
end
if ~strcmp(corr, 'naive')  
    opts.recall = true;
    opts.compute_nulldist = false;
    [FIT_values, FIT_naive] = correction(inputs_1d, outputs, corr, corefunc, opts);
    return
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%             Step 4.B: Compute Probability Distributions                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
entropy_distributions = struct( ...
    'FIT_AB_S', {{'P(S,B_pres,A_past,B_past)'}}, ...
    'FIT_BA_S', {{'P(S,A_pres,B_past,A_past)'}} ...
    );

opts.multidim = true;
required_distributions = {};
for ind = 1:length(indices)
    idx = indices(ind);
    switch possibleOutputs{idx}
        case 'FIT(A->B;S)'
            required_distributions = [required_distributions, entropy_distributions.FIT_AB_S{:}];
        case 'FIT(B->A;S)'
            required_distributions = [required_distributions, entropy_distributions.FIT_BA_S{:}];
    end
end

prob_dists = {};
for i = 1:length(required_distributions)
    switch required_distributions{i}
        case 'P(S,B_pres,A_past,B_past)'
            % B_pres, A_past, B_past
            sourceInput = cat(1, inputs_1d{3},inputs_1d{2},inputs_1d{4});
            prob_dists{i} =  prob_estimator({sourceInput,inputs_1d{end}}, {'P(A,B,C)'}, opts);
        case 'P(S,A_pres,B_past,A_past)'
            % A_pres, B_Past, A_past
            sourceInput = cat(1, inputs_1d{1},inputs_1d{4},inputs_1d{2});
            prob_dists{i} =  prob_estimator({sourceInput,inputs_1d{end}}, {'P(A,B,C)'}, opts);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                  Step 4.C: Compute requested Output Values                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Compute the TE values
FIT_values = cell(1, length(outputs));
atom1 = cell(1, length(outputs));
atom2 = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'FIT(A->B;S)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(S,B_pres,A_past,B_past)')};
            atoms = fit_core(Prob_d, opts);
            FIT_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
        case 'FIT(B->A;S)'
            Prob_d = prob_dists{strcmp(required_distributions, 'P(S,A_pres,B_past,A_past)')};
            atoms = fit_core(Prob_d, opts);
            FIT_values{i} = min(atoms);
            atom1{i} = atoms(1);
            atom2{i} = atoms(2);
    end
end
FIT_naive = FIT_values;
end