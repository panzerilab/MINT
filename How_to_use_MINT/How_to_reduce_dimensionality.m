clc, clear;rng(42);
% -----------------------------------------------------------------------------
% Tutorial: Overview of how to use the MINT toolbox - Dimensionality Reduction
% -----------------------------------------------------------------------------
% In this tutorial, we will explore the concept of dimensionality reduction using a 
% simulated dataset. We will generate a synthetic dataset with a specified number of 
% samples and features, introducing noise. 
% This dataset will be utilized to demonstrate how different techniques, such as 
% Generalized Linear Models (GLM), Support Vector Machines (SVM), and Principal 
% Component Analysis (PCA), can be employed reduce the dimensionality of
% the data and compute information theoretic measures based on the reduced data

nTrials = 600;   
nFeatures = 20;   
noiseLevel = 0.15;  

X = randn(nTrials, nFeatures);   % Random features (normal distribution)
trueBeta = randn(nFeatures, 1); 
zeroIndices = randperm(nFeatures, round(nFeatures * 0.4)); % 40% small beta
trueBeta(zeroIndices) = 0.05;
increaseIndices = randperm(nFeatures, round(nFeatures * 0.3)); % increase 30 % of betas
trueBeta(increaseIndices) = trueBeta(increaseIndices) * 2; 
intercept = 0;                                                % True intercept

% Generate linear combination of features and add intercept
linearCombination = X * trueBeta + intercept;

% Convert linear combination to probabilities using the logistic function
probabilities = 1 ./ (1 + exp(-linearCombination));  % Sigmoid function

% Generate binary labels using these probabilities
labels = binornd(1, probabilities);  % Bernoulli distribution (0 or 1)

% Add noise to the labels (optional, can be turned off)
if noiseLevel > 0
    noisyLabels = labels;
    flipIdx = rand(size(labels)) < noiseLevel;  
    noisyLabels(flipIdx) = 1 - noisyLabels(flipIdx); 
else
    noisyLabels = labels;
end 
% -----------------------------------------------------------------------------
%                               GLM
% -----------------------------------------------------------------------------
% In this section, we utilize the MINT toolbox to perform GLM analysis 
% on the generated data. 

% To fit the GLM model you can use crossvalidation. Therefore ypu can
% specify the field 'cv' as a cell with a String containing the validation
% method and a second cell if needed to specify the number of e.g. Folds
GLM_opts.cv = {'KFold', 3};                       % Cross-validation method: 
                                                  % - 'KFold' (k-fold cross-validation, default k=5)
                                                  % - 'HoldOut' (train-test split), 
                                                  % - 'LeaveMOut' (leave m out),
                                                  % - 'LeaveOut' (leave one out),
                                                  % - 'Resubstitution' (fit on entire dataset).

GLM_opts.glmnet = false;                          % The GLM wrapper supports glmnet. If you want to use it 
                                                  % instead of the matlab native glm function you can set 
                                                  % opts.glmnet to true (default: false)

GLM_opts.distribution = 'binomial';               % Distribution family for GLM:
                                                  % - 'binomial': for binary outcomes (default).
                                                  % - 'normal': for continuous outcomes,
                                                  % - 'poisson': for count data,
                                                  % - 'gaussian': general case for normally distributed responses.
    
GLM_opts.regularization = 'elasticNet';           % Regularization method for GLM:
                                                  % - 'lasso': L1 regularization (alpha = 1)
                                                  % - 'ridge': L2 regularization (alpha = 0.001)
                                                  % - 'elasticNet': combination of L1 and L2 regularization (default).

GLM_opts.weights = [];                            % Sample weights:
                                                  % - Allows for weighted regression where different samples can contribute differently to the loss 
                                                  % function Must be a vector of non-negative values, of the same length as columns of X.
                                                  % (default: [] means everything is weighted with 1).                                               

% If you want to do hyperparamter optimization you have to specify a second
% structure optim_opts that specifies the optimization parameters: 
optim_opts.optimization = true;                   % set to true to perform hyperparameter optimization

optim_opts.optim_cv = {'KFold', 2};               % crossvalidation method for hyperparameter optimization
                                                  % same options as in the outer loop crossvalidation

optim_opts.lambdaRatio = 0.0001;                  % Ratio between the minimum value and maximum value of
                                                  % lambda to generate, if the  parameter "Lambda" is not
                                                  % supplied.  Legal range is [0,1). (default: 0.0001)

optim_opts.NumLambda = 100;                       % Number of lambda values to evaluate:
                                                  % - Specifies how many lambda values should be sampled in the model fitting process.
                                                  % (default: 100).

optim_opts.alphaRange = linspace(0.01, 0.99, 15); % Range of alpha values for Hyperparamter optimization (only used if regularization = elasticNet)
                                                  % - A vector between 0 to 1 (default:linspace(0.01, 0.99, 15))
                                                  % - Alpha controls the balance between L1 and L2 regularization
    
optim_opts.DFmax = Inf;                           % Maximum degrees of freedom:
                                                  % - This parameter limits the complexity of the model.
                                                  % - If set to Inf, it allows any level of model complexity (default: Inf).

optim_opts.RelTol = 1e-4;                         % Relative tolerance for convergence:
                                                  % - A small value that determines the convergence criteria for the optimization algorithm.
                                                  % (default: 1e-4)

optim_opts.MaxIter = 1e4;                          % Maximum number of iterations for the optimization:
                                                  % - This sets the limit on how many iterations the optimization algorithm should run (default: 1e4).

GLM_opts.optim_opts = optim_opts; 

% If you do not want to perform hyperparameter optimization, you can simply provide 
% the values for lambda in the opts field and set optim_opts.optimization = false, 
% or you can omit the optim_opts field entirely in opts. 
% If no optimization is chosen, the alpha value will be set to:
% - 0.01 for Ridge regression,
% - 1 for Lasso regression,
% - 0.9 for ElasticNet.
% Additionally, if no optimization is selected and the lambda field is not provided, 
% the default values for lambda will be [0.001, 0.01, 0.1, 1, 10]
% GLM_opts.lambda = [0.001, 0.01, 0.1, 1, 10];   

outputs = {'all'};                                   % Specify which outputs to return:
                                                     % - 'all' (default) returns all outputs from GLM,
                                                     % - or specify a subset of outputs.
% Function call:                                                     
output_values = glm_wrapper({X', labels},  outputs, GLM_opts);

GLM_labels = output_values{1};                      % Predicted labels.
GLM_labelsBinary = output_values{2}';               % Binary labels.    
GLM_optimized_hyperparameters = output_values{3};   % Optimized hyperparameters (if optimization is performed).            
GLM_betaWeights = output_values{4};                 % Beta weights for the model for each partition.
GLM_intercept_p = output_values{5};                 % Intercept term for the model for each partition.
GLM_meanBetaWeigths = output_values{6};             % Mean beta weights across all partitions.    
GLM_testIdx = output_values{7};                     % Indices of the test set used in each partition.
GLM_confusionMatrix = output_values{8};             % Confusion matrix.  
GLM_diff = abs(GLM_labelsBinary - labels);          % Difference between predicted and true labels.       
GLM_acc = (nTrials-sum(GLM_diff))/nTrials;          % Compute accuracy. 

% Now the predicted labels can be used to compute the population
% information by computing MI of true labels and predicted labels:
MI_opts.bin_method = {'none'};  % labels are already binned
MI_opts.bias = 'shuffSub';      % correct for bias
MI_opts.shuff = 30;
MI_opts.computeNulldist = false;
popInfo_GLM = cell2mat(MI({GLM_labelsBinary', labels'}, {'I(A;B)'}, MI_opts));

% -----------------------------------------------------------------------------
%                               SVM
% -----------------------------------------------------------------------------
% In this section, we utilize the MINT toolbox to perform SVM analysis 
% on the generated data.

% To fit the SVM model, you can use cross-validation. Specify the 'cv' field
% as a cell array with a string containing the validation method and a second
% cell, if needed, to specify the number of folds or splits.
SVM_opts.cv = {'KFold', 3};                       % Cross-validation method:
                                                  % - 'KFold' (k-fold cross-validation, default k=5)
                                                  % - 'HoldOut' (train-test split)
                                                  % - 'LeaveMOut' (leave m out)
                                                  % - 'LeaveOut' (leave one out)
                                                  % - 'Resubstitution' (fit on entire dataset)

SVM_opts.libsvm = false;                          % The SVM wrapper supports libsvm. If you want to use it 
                                                  % instead of MATLAB's native SVM function, set
                                                  % opts.libsvm to true (default: false).

% If you want to perform hyperparameter optimization, you must specify a second
% structure, optim_opts, that contains the optimization parameters.
SVM_optim_opts.optimization = true;               % Set to true to perform hyperparameter optimization.

SVM_optim_opts.optim_cv = {'KFold', 2};           % Cross-validation method for hyperparameter optimization.
                                                  % Same options as in the outer loop cross-validation.

SVM_optim_opts.C_range = [1e-3, 1e3];             % Range of C values for hyperparameter optimization.

SVM_optim_opts.gamma_range = [1e-3, 1e3];         % Range of gamma values for RBF kernel for hyperparameter optimization.

SVM_optim_opts.parallel = true;                   % Whether to run optimization in parallel (default: false).

SVM_optim_opts.MaxIter = 50;                      % Maximum number of iterations for optimization (default: 50).

SVM_opts.optim_opts = SVM_optim_opts;             % Add optimization options to the SVM options.

% If you do not want to perform hyperparameter optimization, you can simply
% provide fixed values for the hyperparameters (opts.hp_C and opts.hp_gamma) and set 
% SVM_optim_opts.optimization = false, or omit the optim_opts field
% entirely. (default C: 1 - default gamma = 1/numFeatures)

outputs = {'all'};                                   % Specify which outputs to return:
                                                     % - 'all' (default) returns all outputs from SVM,
                                                     % - or specify a subset of outputs.
% Function call:
output_values = svm_wrapper({X', labels}, outputs, SVM_opts);

SVM_labels_pred = output_values{1}';              % Predicted labels.
SVM_posteriorProbs = output_values{2};            % Posterior probabilities.
SVM_optimized_hyperparameters = output_values{3}; % Optimized hyperparameters (if optimization is performed).
SVM_betaWeights = output_values{4};               % Beta weights for the model (if linear SVM is used).
SVM_intercept_p = output_values{5};             % Intercept term for the model.
SVM_meanBetaWeigths = output_values{6};           % Mean beta weights across all folds (if linear SVM is used).
SVM_testIdx = output_values{7};                   % Indices of the test set used in each fold.
SVM_diff = abs(SVM_labels_pred - labels);         % Difference between predicted and true labels.
SVM_acc = (nTrials - sum(SVM_diff)) / nTrials;    % Compute accuracy.  

% Now the predicted labels can be used to compute the population
% information by computing MI of true labels and predicted labels:
MI_opts.bin_method = {'none'};  % labels are already binned
MI_opts.bias = 'shuffSub';      % correct for bias
MI_opts.shuff = 30;
popInfo_SVM = cell2mat(MI({SVM_labels_pred', labels'}, {'I(A;B)'}, MI_opts));
%%
% -----------------------------------------------------------------------------
%                               PCA
% -----------------------------------------------------------------------------
% In this section, we utilize the MINT toolbox to perform Principal Component 
% Analysis (PCA) on the generated data.

% To execute PCA, specify the input data as a cell array where the first 
% element contains the feature matrix. The rows represent features, while the 
% columns represent observations.
% For the PCA procedure you can specify additional options:
PCA_opts.standardize = true;                         % Specify whether to standardize the data:
                                                     % - true (default) standardizes the data to have zero mean and unit variance,
                                                     % - false retains the original scale of the data.

% For PCA, you can either specify the number of components to which the data 
% should be reduced or the explained variance, so the data is reduced to the 
% minimum number of components that account for the desired variance. You cannot 
% specify both. If neither is provided, the default is set to an explained 
% variance of 0.95, which will reduce the components to the minimum number 
% needed to explain 95% of the variance:
PCA_opts.numComponents = 4;                          % Set to NaN to select components based on explained variance:
                                                     % - You can specify the number of components directly,
                                                     % - or use explained variance for automatic selection.

PCA_opts.explainedVariance = NaN;                    % Desired amount of variance to explain:
                                                     % - Set between 0 and 1 to indicate the cumulative variance
                                                     %   to be captured by the selected principal components
                                                     %   (default is 0.95).

outputs = {'all'};                                   % Specify which outputs to return:
                                                     % - 'all' (default) returns all outputs from SVM,
                                                     % - or specify a subset of outputs.
% Function call:
output_values = pca_wrapper({X}, outputs, PCA_opts);

% The outputs can include:
% - `coeff`: Principal component coefficients (loadings).
% - `score`: Principal component scores.
% - `latent`: Eigenvalues of the covariance matrix (explained variance).
% - `explained`: Percentage of variance explained by each component.
% - `mean`: Mean of the original data.

PCA_coeff = output_values{1};                       % Coefficients of the principal components.
PCA_score = output_values{2};                       % Scores of the principal components.
PCA_latent = output_values{3};                      % Eigenvalues of the covariance matrix.
PCA_explained = output_values{4};                   % Percentage of variance explained by each component.
PCA_mean = output_values{5};                        % Mean of the original data.

% Now the reduced neural response can be used to compute the population
% information:
MI_opts.bin_method = {'eqpop', 'none'}; 
MI_opts.n_bins = {2};
MI_opts.bias = 'shuffSub';      
MI_opts.shuff = 30;
popInfo_PCA = cell2mat(MI({PCA_coeff', labels'}, {'I(A;B)'}, MI_opts));

% -----------------------------------------------------------------------------
%                               NMF
% -----------------------------------------------------------------------------
% In this section, we utilize the MINT toolbox to perform Non-negative Matrix 
% Factorization (NMF) on the provided data.

% To execute NMF, specify the input data as a cell array where the first 
% element contains the data matrix. The matrix should be N-by-M, where N 
% represents the number of samples and M represents the number of features.
% For the NMF procedure, you can specify additional options:

NMF_opts.newDim = 3;                                  % Desired rank for factorization:
                                                      % - This specifies the target number of reduced dimensions (k).
                                                      % - Must be specified; if not provided, an error will be thrown.

NMF_opts.algorithm = 'als';                           % Algorithm for NMF:
                                                      % - 'als' (default): uses Alternating Least Squares for factorization.
                                                      % - 'mult': uses Multiplicative Update rules.

NMF_opts.w0 = [];                                     % Initial N-by-K matrix for W (optional):
                                                      % - Specify initial values for the factor matrix W.
                                                      % - If left empty, W will be initialized randomly.

NMF_opts.h0 = [];                                     % Initial K-by-M matrix for H (optional):
                                                      % - Specify initial values for the factor matrix H.
                                                      % - If left empty, H will be initialized randomly.

NMF_opts.replicates = 1;                              % Number of times to repeat the factorization:
                                                      % - Determines how many different initializations to run to find the best solution (default is 1).

% Options for further configuration of the algorithm:
options.MaxIter = 100;                                % Maximum number of iterations for convergence:
                                                      % - Sets the limit on iterations for the optimization process (default is 100).
options.Display = 'off';                              % Display options:
                                                      % - 'off': no output (default).
                                                      % - 'final': displays output at the end.
                                                      % - 'iter': displays output at each iteration.
options.TolFun = 1e-4;                                % Tolerance for convergence based on the change in the objective function:
                                                      % - Small values indicate tighter convergence criteria (default is 1e-4).
options.TolX = 1e-4;                                  % Tolerance for convergence based on the change in factor matrices:
                                                      % - Defines the convergence threshold for W and H (default is 1e-4).
options.UseParallel = false;                          % Parallel processing option:
                                                      % - true to enable parallel computation (default is false).
                                                      
options.UseSubstreams = false;                        % Use substreams for reproducibility (default is false).
options.Streams = [];                                 % Set random streams for reproducibility (default is empty).

NMF_opts.options = options;

outputs = {'all'};                                    % Specify which outputs to return:
                                                      % - 'all' (default) returns all outputs from NMF,
                                                      % - or specify a subset of outputs.

% Function call:
output_values = nmf_wrapper({X}, outputs, NMF_opts);

% The outputs can include:
% - `wbest`: Optimal factor matrix W (N-by-K).
% - `hbest`: Optimal factor matrix H (K-by-M).
% - `normbest`: Final residual norm of the factorization.
NMF_wbest = output_values{1};                        % Optimal factor matrix W.
NMF_hbest = output_values{2};                        % Optimal factor matrix H.
NMF_normbest = output_values{3};                     % Final residual norm of the factorization.


% Now the reduced neural response can be used to compute the population
% information:
MI_opts.bin_method = {'eqspace', 'none'};  
MI_opts.n_bins = {2};
MI_opts.bias = 'shuffSub';       
MI_opts.shuff = 30;
popInfo_NMF = cell2mat(MI({NMF_wbest', labels'}, {'I(A;B)'}, MI_opts));
