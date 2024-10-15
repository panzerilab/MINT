function [cFIT_v, shuff_v, cFIT_v_uncorrected] = cFIT(S, X, Y, Z, varargin)
% cFIT calculates the corrected Functional Information Transfer (cFIT) metric.
%
% This function computes the $Z$-conditioned Feature-specific (i.e. $S$-specific)
% Information Transmission (FIT) between a *source* $X$ and a *receiver* $Y$.
%
% Inputs:
% - *S*: must be an array of *nDimsS X nTrials * elements representing the discrete value of the stimulus presented in each trial.
% - *X*: must be an array of *nDimsX X nTimepoints X nTrials* response matrix describing the emitter variables on each of the *nDims* dimensions for each trial.
% - *Y*: must be an array of *nDimsX X nTimepoints X nTrials* response matrix describing the receiver variable on each of the *nDims* dimensions for each trial.
% - *Z*: must be an array of *nDimsX X nTimepoints X nTrials* response matrix describing the the confounding variable on each of the *nDims* dimensions for each trial.
% - *opts*: options used to calculate cFIT (see further notes).
%
% Outputs:
% - *cFIT_biased*: naive, biased estimation of conditional feature-specific information transmission for the input data.
% - *cFIT_unbiased*: unbiased estimate of feature-specific information transmission.  This is returned as second argument only if `opts.bias` is not `'naive'`.
% - *cFIT_pval*: p-value of the measure tested against the null hypothesis of Y being conditionally independent from X given S
%
% Further notes:
% The *opts* structure can have the following fields:
%
% | field                                | description                                                                                                                                                                                                                                             | allowed values                                                                                                                                                                                                                                                               | default   |
% |--------------------------------------|---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------|-----------|
% | opts.bias                            | specifies the bias correction method                                                                                                                                                                                                                    | `'naive'` (no bias correction)<br>`'le'` (linear extrapolation)<br>`'qe'` (quadratic extrapolation)                                                                                                                                                                          | `'naive'` |
% | opts.bin_methodX                     | specifies the binning method for the X (emitter) signals                                                                                                                                                                                                | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details    | `'eqpop'` |
% | opts.bin_methodY                     | specifies the binning method for the Y (receiver) signals                                                                                                                                                                                               | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details    | `'eqpop'` |
% | opts.bin_methodZ                     | specifies the binning method for the Z (confounding) signals                                                                                                                                                                                            | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details    | `'eqpop'` |
% | opts.bin_methodS                     | specifies the binning method for the stimulus signals                                                                                                                                                                                                   | `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details    | `'none'`  |
% | opts.n_binsX                         | number of bins to be used to reduce the dimensionality of the response X                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                                                      | 3         |
% | opts.n_binsY                         | number of bins to be used to reduce the dimensionality of the response Y                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                                                      | 3         |
% | opts.n_binsZ                         | number of bins to be used to reduce the dimensionality of the response Y                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                                                      | 3         |
% | opts.n_binsS                         | number of bins to be used to reduce the dimensionality of the stimulus S                                                                                                                                                                                | int > 1                                                                                                                                                                                                                                                                      | 2         |
% | opts.nullhyp                         | If 'true' test the cFIT value against the null hypothesis of Y being conditionally independent from X once S is given, via a permutation test. Otherwise, don't test significance                                                                        | bool                                                                                                                                                                                                                                                                         | false     |
% | opts.nh_perm                         | number of permutations used to build the null hypothesis distribution, by default 200
% | opts.taux                            | lag considered for the causing signal. Lag is specified as **strictly negative** integer.                                                                                                                                                               | int < 0                                                                                                                                                                                                                                                         | 3         |
% | opts.tauy                            | lag considered for the caused signal. Lag is specified as **strictly negative** integer. 
% | opts.tauz                            | lag considered for the confounding signal. Lag is specified as **strictly negative** integer.         
% | opts.tpres                           | timepoint which should be considered as present in the computation (integer)


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
        [outputs, opts] = check_inputs('cFIT',inputs,varargin{:});
    end
else
    [outputs, opts] = check_inputs('cFIT',inputs,varargin{:});
end

%% TO BE DONE

% set default options
if nargin == 4
    opts = default_opts;
else
    opts = varargin{1};
    missing_fields = setdiff(fieldnames(default_opts), fieldnames(opts));
    for i=1:size(missing_fields)
        opts.(missing_fields{i}) = default_opts.(missing_fields{i});
    end
end

if ~any(arrayfun(@(x) strcmp(x, opts.bias), ["qe", "naive"]))                           %tbd: adding more bias
    msg = "User defined bias correction method: `" + opts.bias...
        + "` cannot be used for cFIT.";
    error('cFIT:UndefinedBiasCorrectionMethod', msg);
end

allowedMethods = ["none", "eqpop", "eqspace", "ceqspace", "gseqspace"];
binMethods = {opts.bin_methodX, opts.bin_methodY, opts.bin_methodZ};
methodNames = {'bin_methodX', 'bin_methodY', 'bin_methodZ'};
for i = 1:length(binMethods)
    currentMethod = binMethods{i};
    currentName = methodNames{i}; 
    if ~any(arrayfun(@(x) strcmp(x, currentMethod), allowedMethods))
        if ~exist(currentMethod, 'file')
            msg = "User defined binning method: `" + currentMethod + "` for " + currentName + " cannot be found in the current MATLAB search path.";
            error('cFIT:InvalidBinningMethod', msg);
        elseif exist(currentMethod, 'file') ~= 2
            msg = "User defined binning method: `" + currentMethod + "` for " + currentName + " could be found in the current MATLAB search path but does not appear to be a MATLAB function.";
            error('cFIT:InvalidBinningMethod', msg);    
        end
    end
end


opts.redundancy_measure = 'I_min';
% opts = check_opts(opts, 'cFIT');
opts.n_trials =n_trials;

% bin neural responses if required
if ~strcmp(opts.bin_methodS, 'none')
    binr_opts.nb = opts.n_binsS;
    binr_opts.method = opts.bin_methodS;
    binr_opts.supressWarnings = opts.supressWarnings;
    S_b = binr(S, binr_opts) + 1;
else
    S_b = S - min(S) + 1;
end


xtau = opts.taux;
ytau = opts.tauy;
ztau = opts.tauz;

tP_present = opts.tpres;
tP_past_X = tP_present - abs(xtau);
tP_past_Y = tP_present - abs(ytau);
tP_past_Z = tP_present - abs(ztau);

Y_all = Y; 

hX = squeeze(X(:, tP_past_X, :));
hY = squeeze(Y(:, tP_past_Y, :));
hZ = squeeze(Z(:, tP_past_Z, :));
Y = squeeze(Y(:, tP_present, :));

hX = reshape(hX, [size(X,1), size(X,3)]);
hY = reshape(hY, [size(Y_all,1), size(Y_all,3)]);
Y = reshape(Y, [size(Y_all,1), size(Y_all,3)]);
hZ = reshape(hZ, [size(Z,1), size(Z,3)]);

if ~strcmp(opts.bin_methodX, 'none')
    binr_opts.nb = opts.n_binsX;
    binr_opts.method = opts.bin_methodX;
    binr_opts.supressWarnings = opts.supressWarnings;
    hX_b = binr(hX, binr_opts) + 1; % The '+ 1' is added because zeroes in the array would cause an error in the probabilityDist function
else
    hX_b = hX - min(hX) + 1;
end

if ~strcmp(opts.bin_methodY, 'none')
    binr_opts.nb = opts.n_binsY;
    binr_opts.method = opts.bin_methodY;
    binr_opts.supressWarnings = opts.supressWarnings;
    hY_b = binr(hY, binr_opts) + 1;
    Y_b = binr(Y, binr_opts) + 1;
else
    hY_b = hY - min(hY) + 1;
    Y_b = Y - min(Y) + 1;
end

if ~strcmp(opts.bin_methodZ, 'none')
    binr_opts.nb = opts.n_binsZ;
    binr_opts.method = opts.bin_methodZ;
    binr_opts.supressWarnings = opts.supressWarnings;
    hZ_b = binr(hZ, binr_opts) + 1;
else
    hZ_b = hZ - min(hZ) + 1;
end

% if multi-dimensional response map it onto a 1D space
if length(hX_b(:,1)) > 1
    hX_b =reduce_dim(hX_b, 1);
end
if length(hY_b(:,1)) > 1
    hY_b = reduce_dim(hY_b, 1);
    Y_b = reduce_dim(Y_b,1);
end
if length(S_b(:,1)) > 1
    S_b = reduce_dim(S_b,1);
end
if length(hZ_b(:,1)) > 1
    hZ_b = reduce_dim(hZ_b,1);
end
opts.method='none';
cFIT_v_uncorrected = 0;
if strcmp(opts.bias,'naive')
    [p_S4] = pdf([Y_b; hX_b; hY_b; hZ_b; S_b], opts);
    cFIT_v = cfit_core(p_S4, opts);
else
    [cFIT_v, cFIT_v_uncorrected] = extrapolation(@cfit_core, [Y_b; hX_b; hY_b; hZ_b; S_b], opts, false);
end


shuff_v = zeros(length(opts.shuff_variables),opts.shuff);
for s = 1:length(opts.shuff_variables)
    % initialize arrays for shuff results
    switch opts.shuff_variables{s}
        case 'S'
            shuffled_var = S_b;
        case 'hX'
            shuffled_var = hX_b;
        case 'hY'
            shuffled_var = hY_b;
        case 'Y'
            shuffled_var = Y_b;
    end
    switch opts.conditioning_vars{s}
        case 'all'
            % un-conditioned case
            cond_var = [];
            consistent = 0;
        case 'Sconditioned'
            % S conditioned
            cond_var = S_b;
            consistent = 1;
        case 'hXconditioned'
            % hX conditioned
            cond_var = hX_b;
            consistent = 1;
        case 'hYconditioned'
            % hY conditioned
            cond_var = hY_b;
            consistent = 1;
        case 'Yconditioned'
            % Y conditioned
            cond_var = Y_b;
            consistent = 1;
    end

    for b = 1:opts.shuff
        shuffled_var_tmp = shuffle(cond_var', shuffled_var', consistent, 1);
        shuffled_var = shuffled_var_tmp';

        switch opts.shuff_variables{s}
            case 'S'
                [p_S4] = pdf([Y_b; hX_b; hY_b; hZ_b; S_b], opts);
            case 'hX'
                [p_S4] = pdf([Y_b; hX_b; hY_b; hZ_b; S_b], opts);
            case 'hY'
                [p_S4] = pdf([Y_b; hX_b; hY_b; hZ_b; S_b], opts);
            case 'Y'
                [p_S4] = pdf([Y_b; hX_b; hY_b; hZ_b; S_b], opts);
        end
        shuff_v(s,d) = cfit_core(p_S4, opts);
    end
end


end