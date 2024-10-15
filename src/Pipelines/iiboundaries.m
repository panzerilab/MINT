function [II_Info, theta, bw_stim, bw_choice] = iiboundaries(S, R, C, nfolds, svm_opts, ii_opts)
%%% *function [II_Info, theta] = iiboundaries(S, R, C, nfolds, svm_opts, ii_opts)*
%%%
%%% ### Description
%%%  This function computes the intersection information (II_Info) and the angle
%%%  (theta) between the decision boundaries of stimuli and choices based on neural
%%%  responses.
%%%
%%% ### Inputs:
%%%   - S: Binary vector of stimulus conditions (1 x nTrials)
%%%   - R: Matrix of neural responses (nNeurons x nTrials)
%%%   - C: Binary vector of choice conditions (1 x nTrials)
%%%   - nfolds: Number of cross-validation folds
%%%   - svm_opts: Structure with optional parameters for configuring the SVM model and optimization process
%%%       - optimize_params: Flag for hyperparameter optimization (default is false)
%%%       - cv_type: Cross-validation type, 'KFold' or 'LeaveOneOut'
%%%       - K: Number of folds for KFold cross-validation (required if cv_type is 'KFold')
%%%       - optim_reps: Number of repetitions for hyperparameter optimization (default is 10)
%%%       - svm_family: Type of SVM kernel, 'linear' or 'RBF' for Radial Basis Function
%%%       - hp_C: Hyperparameter C for SVM regularization (default is 1)
%%%       - hp_gamma: Hyperparameter gamma for RBF kernel (default is 1/n_features)
%%%       - libsvm: Flag for using LIBSVM (true) or FITCSVM (false) (default is false)
%%%   - ii_opts: Structure with optional parameters for intersection information
%%%       - bias: specifies the bias correction method   
%%%       - max_draws_per_split_number: specifies the maximum number of draws on which to calculate the unbiased estimate of II. E.g. `opts.max_draws_per_split_number = 10`, the II is calculated as average on 20 trials for 2 splits and 10 trials for 1 split. (ignored if opts.bias = 'naive') 
%%%       - bin_methodS/opts.bin_methodR/opts.bin_methodC: `'none'` (no binning)<br>`'eqpop'` (evenly populated binning)<br>`'eqspace'`(evenly spaced binning)<br>`'ceqspace'`(centered evenly spaced binning)<br>`'gseqspace'`(gaussian evenly spaced binning)<br> See the documentation of [[binr|binr]] function for more details
%%%       - n_binsS/opts.n_binsR/opts.n_binsC: number of bins to be used to reduce the dimensionality
%%%       - btsp: number of bootstrap operations to perform for significance testing (those will be performed independently on each of the variables listed in `opts.btsp_variables`) 
%%%       - btsp_variables: list of variables to be bootstrapped, specified as a cell array of strings. If multiple bootstrapping operations are requested, `opts.btsp_variables` can contain multiple variables. In this case a `opts.btsp_type` option should be specified for each `opts.btsp_variables`
%%%       - btsp_type: type of bootstrapping to be applied for each variable (`'all'` shuffles all values of the corresponding variable in `opts.btsp_variables` across trials, while `'$VAR$conditioned'` shuffles trials by conditioning on values of the variable specified in the substring `$VAR$`)|
%%%       - btsp_bias: bias correction method to be applied for each variable (same values than in opts.bias)    
%%%
%%% ### Outputs:
%%%   - II_Info; Intersection information between stimulus, neural responses and choice
%%%   - theta - Angle between decision boundaries of stimulus and choice


sanity_check(R);
sanity_check(S);
sanity_check(C);

nTrials = length(S);
nNeurons = size(R, 1);

cvPartition = cvpartition(S, 'KFold', nfolds); 
bW_stim = zeros(nfolds, nNeurons);
ic_stim = zeros(nfolds, 1);

bW_choice = zeros(nfolds, nNeurons);
ic_choice = zeros(nfolds, 1);

Stim_dec = zeros(1, nTrials);
Choice_dec = zeros(1, nTrials);

for cvRep = 1:nfolds
    [Stim_dec(1,cvPartition.test(cvRep)) , ~, ~, ~, ~,bW_stim(cvRep,:), ic_stim(cvRep)] = svm_pipeline(R',S', find(cvPartition.test(cvRep)), svm_opts); 
    [Choice_dec(1,cvPartition.test(cvRep)), ~, ~, ~, ~,bW_choice(cvRep,:), ic_choice(cvRep)] = svm_pipeline(R',C',find(cvPartition.test(cvRep)), svm_opts); 
end 

dec_SC=map_Nd_array_to_1d([Stim_dec; Choice_dec]);
II_Info = II(Stim_dec,dec_SC, Choice_dec,ii_opts); 

if nNeurons == 2
    slope_s = mean(-bW_stim (:,1)./bW_stim(:,2));
    slope_c = mean(-bW_choice(:,1)./bW_choice(:,2));
    theta = atan(abs((slope_s - slope_c) / (1 + slope_c * slope_s)));
    theta = rad2deg(theta);
else
    normal_vector_s = mean(bW_stim);   
    normal_vector_c = mean(bW_choice);
    cos_theta = dot(normal_vector_s, normal_vector_c) / (norm(normal_vector_s) * norm(normal_vector_c));
    theta = acos(cos_theta);
    theta = rad2deg(theta);
end
end