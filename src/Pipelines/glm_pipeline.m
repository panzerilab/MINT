function [predictedLabels,posteriorProbs, labelsTest, dataTest, confMat,  varargout] = glm_pipeline(data, labels, testSet, opts)                                               
% glm_pipeline: Train and evaluate Generalized Linear Model (GLM) using lasso,ridge or elastic net regularization and either LASSOGLM or (CV)GLMNET.
%   This function allows for the training and testing of GLM models with optional regularization methods, including lasso, ridge, or elastic net. 
%   It provides flexibility in choosing regularization types and additional options for k-fold lambda cross-validation 
%
% Input:
%   - data: Training data matrix where each row represents a sample, and each column represents a feature.
%   - labels: True labels corresponding to the training data samples.
%   - testSet: Indices of samples to be used as the test set.
%   - opts: Structure containing optional parameters for configuring the GLM model.
%       - regularization: Type of regularization to be applied ('lasso', 'ridge', or 'elasticNet').
%       - alpha: Elastic net mixing parameter (used if regularization is 'elasticNet', default is 0.95).
%       - CV: Flag indicating whether to perform cross-validation for regularization parameter selection (default is false).
%       - folds: Number of folds for cross-validation (required if CV is true).
%       - lambda: Regularization parameter values (used if CV is false, default is [0.01, 0.1, 1, 10]).
%       - glmnet: Flag indicating whether to use the glmnet implementation for regularization (default is false).
%
% Output:
%   - predictedLabels: Predicted labels for the test set.
%   - labelsTest: True labels for the test set.

%% Check input parameters

% Validate input parameters and set default values if necessary
if any(isnan(data(:))) || any(isnan(labels)) || any(isnan(testSet(:)))
    error("One or more variables contain NaNs. Aborting.")
end

if ~isfield(opts,'distribution')
    opts.distribution = 'normal';
end

if strcmp(opts.regularization, 'lasso')
    opts.alpha = 1;
end 

if strcmp(opts.regularization, 'ridge')
    opts.alpha = 0.001;
end 

if strcmp(opts.regularization, 'elasticNet')
    if ~isfield(opts,'alpha')
        opts.alpha = 0.95;
    end
    if ~isfield(opts,'CV')
        opts.CV = false;
    end
end 

if ~isfield(opts,'glmnet')
    opts.glmnet = false;
end

if opts.CV == true
    if ~isfield(opts,'folds')
        opts.folds = 3;
    end
end

if opts.CV == false
    if ~isfield(opts,'lambda')
        opts.lambda = [0.001, 0.01, 0.1, 1, 10, 100];
    end
end

%% Select training and test set

% Divide the data into training and test sets based on user-defined testSet
nTrials = length(labels);
nTrialsTest = length(testSet);
nTrialsTrain = nTrials-nTrialsTest;

if nTrialsTrain == 0 || nTrialsTest == 0
    error('The number of trials for test or train cannot be zero. Modify your test fraction.');
end

train_test_partition.test = testSet;
train_test_partition.training = setdiff(1:nTrials, testSet);

% Training and test set for main computation
dataTrain = data(train_test_partition.training',:);
labelsTrain = labels(train_test_partition.training');
dataTest = data(train_test_partition.test',:);
labelsTest = labels(train_test_partition.test');


if opts.CV == true
    if opts.glmnet == true
        cvfit = cvglmnet(dataTrain, labelsTrain, opts.distribution, struct('alpha', opts.alpha), [], opts.folds);
        varargout{1} = cvfit.lambda_min;
        predictedLabels = cvglmnetPredict(cvfit, dataTest, cvfit.lambda_min, 'response');
    else
        [B, FitInfo] = lassoglm(dataTrain, labelsTrain, opts.distribution, 'Alpha', opts.alpha, 'CV', opts.folds);
        B0 = FitInfo.Intercept(FitInfo.IndexMinDeviance);
        coef = [B0; B(:, FitInfo.IndexMinDeviance)];
        varargout{1} = FitInfo.LambdaMinDeviance;
        predictedLabels = glmval(coef, dataTest, 'logit');
        %predictedLabels_binom = (predictedLabels>=0.5);
        %confMat = confusionchart(labelsTest,double(predictedLabels_binom));
    end
else
    if opts.glmnet == true
        fit = glmnet(dataTrain, labelsTrain, opts.distribution, struct('alpha', opts.alpha, 'lambda', opts.lambda));
        [~, idxBestLambda] = max(fit.dev);
        varargout{1} = fit.lambda(idxBestLambda);
        predictedLabels = glmnetPredict(fit, dataTest, fit.lambda(idxBestLambda), 'response');
    else
        [B, FitInfo] = lassoglm(dataTrain, labelsTrain, opts.distribution, 'Alpha', opts.alpha, 'Lambda', opts.lambda);
        [~, idxBestLambda] = min(FitInfo.Deviance);
        varargout{1} = FitInfo.Lambda(idxBestLambda);   
        B0 = FitInfo.Intercept(idxBestLambda);
        coef = [B0; B(:,idxBestLambda)];
        predictedLabels = glmval(coef,dataTest,'logit');
        %predictedLabels_binom = (predictedLabels>=0.5);
        %confMat = confusionchart(labelsTest,double(predictedLabels_binom));
    end
end

end

