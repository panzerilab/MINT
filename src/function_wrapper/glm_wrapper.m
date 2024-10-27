function output_values = glm_wrapper(input,varargin)
% *function output_values = glm_wrapper(input,outputs, opts)*
% GLM_pipeline - Generalized Linear Model Implementation
%
% This function implements a Generalized Linear Model (GLM)
% classification, offering options for regularization techniques such as Elastic Net,
% as well as cross-validation and hyperparameter optimization.
%
% Usage:
%   output_values = GLM_pipeline(input, varargin)
%
% Inputs:
%   - input: A cell array containing:
%       - input{1}: Feature matrix (data), where rows represent observations and columns
%         represent features.
%       - input{2}: Label vector (labels) corresponding to the observations in the feature
%         matrix.
%
%   - varargin: Optional arguments, which may include:
%       - outputs: A cell array specifying which outputs to return (default is {'all'}).
%       - opts: A structure containing options, including:
%           - glmnet: Boolean indicating whether to use glmnet for regularization (default false).
%           - distribution: Distribution family for the GLM (default is 'normal').
%             Supported options: 'normal', 'binomial', 'poisson', 'gamma', 'inverse gaussian'.
%           - regularization: Regularization method to use ('elasticNet' is default).
%             Supported options: 'lasso', 'ridge', 'elasticNet'.
%           - alpha: Elastic Net mixing parameter (default is 0.9 for elasticNet).
%             This balances between ridge (alpha=0) and lasso (alpha=1).
%           - lambda: Regularization strength, can be a scalar or vector of values for hyperparameter search (Default: [0.001, 0.01, 0.1, 1, 10]).
%           - cv_type: Cross-validation type (default is 'KFold'). Supported options: 'KFold', 'HoldOut', 'LeaveMOut', 'LeaveOneOut', 'Resubstitution'.
%           - K_folds: Number of folds for K-fold cross-validation (default is 5).
%           - weights: Optional observation weights (default is all weighted with 1).
%           - CV: Boolean indicating whether to perform cross-validation for model evaluation (default: false).
%           - optim_opts: Structure for hyperparameter optimization options, including:
%               - optimize: Boolean indicating whether to enable hyperparameter optimization (default false).
%               - lambdaRange: Range of lambda values for hyperparameter search (default is [0.001, 0.01, 0.1, 1, 10, 100]).
%               - alphaRange: Vector specifying range of alpha values to search (default is linspace(0.01, 0.99, 15)).
%               - NumLambda: Number of lambda values to evaluate in hyperparameter search (default is 100).
%               - lambdaRatio: Ratio of the largest to smallest lambda in the search (default is 0.0001).
%               - DFmax: Maximum number of non-zero coefficients allowed (default is Inf).
%               - RelTol: Relative tolerance for model fitting (default is 1e-4).
%               - MaxIter: Maximum number of iterations for fitting (default is 1e4).
%               - optim_cv: Cell array specifying cross-validation method for optimization (default is {'Resubstitution'}).
%                 Supported options: 'KFold', 'HoldOut', 'LeaveMOut', 'LeaveOneOut', 'Resubstitution'.
%
% Outputs:
%   - output_values: A cell array containing the requested output values, which may include:
%       - 'labels': Predicted class labels for the test data.
%       - 'labelsBinary': Binary predictions (for binomial models) based on a threshold of 0.5.
%       - 'optimized_hyperparameters': Optimal lambda and alpha values selected via cross-validation.
%       - 'confusionMatrix': Confusion matrix of predictions for each fold.
%       - 'betaWeights': Coefficients of the GLM model for each fold.
%       - 'meanBetaWeights': Mean of the coefficients across all cross-validation folds.
%       - 'testIdx': Indices of the test data in each fold.
%       - 'intercept': Intercept terms for the model in each fold.
%
% Note: Hyperparameter optimization is supported via cross-validation, allowing the model 
% to tune parameters such as lambda and alpha. Users should configure `optim_opts` to enable
% this functionality. The `opts` structure allows for extensive customization of the model fitting.
%
% Example:
%   input = {dataMatrix, labelVector};
%   opts = struct('glmnet', false, 'distribution', 'binomial', 'cv_type', 'KFold', 'K_folds', 5);
%   outputs = {'labels', 'confusionMatrix'};
%   output_values = glm_wrapper(input, outputs, opts);

% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT.
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses>.



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check input, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultOpts.glmnet = false;
defaultOpts.distribution = 'normal';
defaultOpts.regularization = 'elasticNet';
defaultOpts.lambda = [0.001, 0.01, 0.1, 1, 10];
defaultOpts.cv_type = {'KFold', 5};
defaultOpts.weights = [];

defaultOptimOpts.lambdaRatio = 0.0001;
defaultOptimOpts.alphaRange = linspace(0.01, 0.99, 15);
defaultOptimOpts.NumLambda = 100;
defaultOptimOpts.DFmax = Inf;
defaultOptimOpts.RelTol = 1e-4;
defaultOptimOpts.MaxIter = 1e4;



if nargin < 1
    msg = 'Not enough input arguments. See `help glm_wrapper` for usage info';
    error('GLM:notEnoughInput', msg);
elseif nargin == 1
    opts = defaultOpts;
    outputs = {'all'};
    opts.isChecked = true;
elseif nargin == 2
    if iscell(varargin{1})
        outputs = varargin{1};
        opts = defaultOpts;
        opts.isChecked = true;
    elseif isstruct(varargin{1})
        opts = varargin{1};
        outputs = {'all'};
        opts.isChecked = false;
    end
elseif nargin == 3
    opts = varargin{2};
    outputs = varargin{1};
    opts.isChecked = false;
end

if ~isfield(opts, 'NaN_handling')
    opts.NaN_handling = 'error';
end
input = nan_method(input, opts.NaN_handling);

data = input{1};
labels = input{2};

if ~opts.isChecked
    default_fields= fieldnames(defaultOpts);
    is_field_present = ismember(default_fields, fieldnames(opts));
    missing_fields = default_fields(~is_field_present);
    for i=1:size(missing_fields,1)
        missing_field_name = missing_fields{i};
        opts.(missing_fields{i}) = defaultOpts.(missing_fields{i});
    end
end

if ~isfield(opts,'optim_opts')
    optimization_opts.optimize = false;
else
    optimization_opts = opts.optim_opts;
    if ~isfield( optimization_opts, 'optimize')
        optimization_opts.optimize = true;
    end
    default_fields= fieldnames(defaultOptimOpts);
    is_field_present = ismember(default_fields, fieldnames(optimization_opts));
    missing_fields = default_fields(~is_field_present);
    for i=1:size(missing_fields,1)
        missing_field_name = missing_fields{i};
        optimization_opts.(missing_fields{i}) = defaultOptimOpts.(missing_fields{i});
    end
end
if strcmp(opts.regularization, 'lasso')
    opts.alpha = 1;
elseif strcmp(opts.regularization, 'ridge')
    opts.alpha = 0.001;
elseif strcmp(opts.regularization, 'elasticNet') && ~optimization_opts.optimize
    if ~isfield(opts,'alpha')
        opts.alpha = 0.9;
    end
end

if optimization_opts.optimize
    if strcmp(opts.regularization, 'elasticNet')
        if ~isfield(opts,'alpha')
            opts.alpha = optimization_opts.alphaRange;
        end
    end
    if ~isfield(optimization_opts,'optim_cv')
        optimization_opts.optim_cv = {'Resubstitution'};
        warning('No cross-validation method specified for the hyperparameter optimization. Using default method: ''Resubstitution''. The GLM will use the entire dataset to optimize the hyperparameters.');
    end
    optim_cv_method = string(optimization_opts.optim_cv{1});
    valid_methods = ["KFold", "HoldOut", "LeaveMOut", "LeaveOut", "Resubstitution"];
    if ~ismember(optim_cv_method, valid_methods)
        error('GLM:The chosen cross-validation method for hyperparameter optimization is not acceptable');
    end
    if strcmp(optim_cv_method,'KFold')
        if length(optimization_opts.optim_cv) > 1
            optim_cv_partition =  optimization_opts.optim_cv{2};
        else
            optim_cv_partition =  3;
            warning('Using default partition value for crossvalidation: %d.', optim_cv_partition);
        end
    elseif strcmp(optim_cv_method,'HoldOut')
        if length(optimization_opts.optim_cv) > 1
            optim_cv_partition =  optimization_opts.optim_cv{2};
        else
            optim_cv_partition =  0.2;
            warning('Using default partition value for crossvalidation: %d.', optim_cv_partition);
        end
    elseif strcmp(optim_cv_method,'LeaveMOut')
        if length(optimization_opts.optim_cv) > 1
            optim_cv_partition =  optimization_opts.optim_cv{2};
        else
            optim_cv_partition =  2;
            warning('Using default partition value for crossvalidation: %d.', optim_cv_partition);
        end
    else
        optim_cv_partition =  NaN;
    end
end

if ~isfield(opts,'cv')
    opts.cv = {'Resubstitution'};
    warning('No cross-validation method specified. Using default method: ''Resubstitution''. The GLM will be fitted on the entire dataset without cross-validation.');
end

if iscell(opts.cv)
    cv_method =  opts.cv{1};
    valid_methods = ["KFold", "HoldOut", "LeaveMOut", "LeaveOut", "Resubstitution"];
    if ~ismember(cv_method, valid_methods)
        error('GLM:The chosen cross-validation method for hyperparameter optimization is not acceptable');
    end
    if strcmp(cv_method,'KFold')
        if length(opts.cv) > 1
            cv_partition =  opts.cv{2};
        else
            cv_partition =  5;
            warning('Using default partition value for crossvalidation: %d.', cv_partition);
        end
        cvPartition = cvpartition(input{2}, cv_method, cv_partition);
    elseif strcmp(cv_method,'HoldOut')
        if length(opts.cv) > 1
            cv_partition =  opts.cv{2};
        else
            cv_partition =  0.2;
            warning('Using default partition value for crossvalidation: %d.', cv_partition);
        end
        cvPartition = cvpartition(input{2}, cv_method, cv_partition);
    elseif strcmp(cv_method,'LeaveMOut')
        if length(opts.cv) > 1
            cv_partition =  opts.cv{2};
        else
            cv_partition =  2;
            warning('Using default partition value for crossvalidation: %d.', cv_partition);
        end
        cvPartition = cvpartition(input{2}, cv_method, cv_partition);
    else
        cvPartition = cvpartition(input{2}, cv_method);
    end

else
    error('GLM:Invalid input for cross-validation. The cross-validation option must be provided as a cell array. For ''KFold'', specify the number of folds as {''KFold'', numFolds}, e.g., {''KFold'', 5}. The other valid option is {''LeaveOneOut''}.');
end


possibleOutputs = {'labels', 'labelsBinary', 'optimized_hyperparameters', 'betaWeights', 'intercept','meanBetaWeigths', 'testIdx','confusionMatrix'};
if ismember('all', outputs)
    outputs = possibleOutputs;
end
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('GLM:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                       Step 3: Select training and test set                    %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predictedLabels = zeros(length(labels), 1);
predictedLabels_binom = zeros(length(labels), 1);
confMat.partition = cell(cvPartition.NumTestSets, 1);
partition_betaWeights = zeros(size(input{1}, 1),cvPartition.NumTestSets);
optimHyp = struct();
intercept = zeros(1, cvPartition.NumTestSets);
testIdx = cell(cvPartition.NumTestSets, 1);
for partition = 1:cvPartition.NumTestSets
    test_idx = find(cvPartition.test(partition));
    nTrials = length(labels);
    nTrialsTest = length(test_idx);
    nTrialsTrain = nTrials-nTrialsTest;
    if nTrialsTrain == 0 || nTrialsTest == 0
        train_test_partition.test = test_idx;
        train_test_partition.training = test_idx;
    else
        train_test_partition.test = test_idx;
        train_test_partition.training = setdiff(1:nTrials, test_idx);
    end
    testIdx{partition} = train_test_partition.test;
    dataTrain = data(:,train_test_partition.training)';
    labelsTrain = labels(train_test_partition.training);
    dataTest = data(:,train_test_partition.test)';
    labelsTest =labels(train_test_partition.test);
    if optimization_opts.optimize == true
        if isnan(cv_partition)
            optim_cvPartition =  cvpartition(labelsTrain, optim_cv_method);
        else
            optim_cvPartition =  cvpartition(labelsTrain, optim_cv_method, cv_partition);
        end
        if opts.glmnet == true
            cvfit = cvglmnet(dataTrain, labelsTrain, opts.distribution, struct('alpha', opts.alpha), [], opts.CV_folds);
            optimHyp.partition = cvfit.lambda_min;
            predictedLabels(test_idx) = cvglmnetPredict(cvfit, dataTest, cvfit.lambda_min, 'response');
        else
            if length(opts.alpha) > 1
                bestDeviance = Inf;
                for i = 1:length(opts.alpha)
                    alpha = opts.alpha(i);
                    [B, FitInfo] = lassoglm(dataTrain, labelsTrain, opts.distribution, ...
                        'Alpha', alpha, 'LambdaRatio', optimization_opts.lambdaRatio, 'NumLambda', optimization_opts.NumLambda, ...
                        'CV', optim_cvPartition, 'DFmax', optimization_opts.DFmax, 'RelTol', optimization_opts.RelTol, 'MaxIter', optimization_opts.MaxIter, 'Weights', opts.weights);
                    [~, idxMinDeviance] = min(FitInfo.Deviance);
                    if FitInfo.Deviance(idxMinDeviance) < bestDeviance
                        bestDeviance = FitInfo.Deviance(idxMinDeviance);
                        bestAlpha = alpha;
                        bestLambda = FitInfo.Lambda(idxMinDeviance);
                        bestCoef = [B(:, idxMinDeviance)];
                        bestIntercept = FitInfo.Intercept(idxMinDeviance);
                    end
                end
                optimHyp.lambda(partition) = bestLambda;
                optimHyp.alpha(partition) = bestAlpha;
                partition_betaWeights(:, partition) = bestCoef;
                intercept(partition) = bestIntercept;
            else
                alpha = opts.alpha;
                [B, FitInfo] = lassoglm(dataTrain, labelsTrain, opts.distribution, ...
                    'Alpha', alpha, 'LambdaRatio', optimization_opts.lambdaRatio, 'NumLambda', optimization_opts.NumLambda, ...
                    'CV', optim_cvPartition, 'DFmax', optimization_opts.DFmax, 'RelTol', optimization_opts.RelTol, 'MaxIter', optimization_opts.MaxIter,'Weights', opts.weights);
                [~, idxMinDeviance] = min(FitInfo.Deviance);
                optimHyp.partition = FitInfo.Lambda(idxMinDeviance);
                bestCoef = [B(:, idxMinDeviance)];
                intercept(partition) = FitInfo.Intercept(idxMinDeviance);
                partition_betaWeights(:, partition) = bestCoef;
            end
            predictedLabels(test_idx) = glmval([FitInfo.Intercept(idxMinDeviance); bestCoef], dataTest, 'logit');
            predictedLabels_binom(test_idx) = (predictedLabels(test_idx) >= 0.5);
            if ismember('confusionMatrix', outputs)
                fig = figure('Visible', 'off'); 
                confMat.partition{partition} = confusionchart(labelsTest, double(predictedLabels_binom(test_idx)), 'Visible', 'off');
                delete(fig);
            end 
        end
    else
        alpha = opts.alpha;
        lambda = opts.lambda;
        [B, FitInfo] = lassoglm(dataTrain, labelsTrain, opts.distribution, 'Alpha', alpha, 'Lambda', lambda, 'Weights', opts.weights);
        if length(lambda) > 1
            [~, idxBestLambda] = min(FitInfo.Deviance);
            intercept(partition) = FitInfo.Intercept(idxBestLambda);
            coef = [B(:, idxBestLambda)];
            optimHyp.partition = FitInfo.Lambda(idxBestLambda);
        else
            intercept(partition) = FitInfo.Intercept;
            coef = B;
            optimHyp.partition = lambda;
        end
        predictedLabels(test_idx) = glmval([intercept(partition); coef], dataTest, 'logit');
        predictedLabels_binom(test_idx) = (predictedLabels(test_idx) >= 0.5);
        if ismember('confusionMatrix', outputs)
            fig = figure('Visible', 'off');
            confMat.partition{partition} = confusionchart(labelsTest, double(predictedLabels_binom(test_idx)), 'Visible', 'off');
            delete(fig);
        end
    end
end

output_values = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'labels'
            output_values{i} = predictedLabels';
        case 'labelsBinary'
            output_values{i} = predictedLabels_binom';
        case 'optimized_hyperparameters'
            output_values{i} = optimHyp;
        case 'confusionMatrix'
            output_values{i} =  confMat;
        case 'betaWeights'
            output_values{i} = partition_betaWeights;
        case 'meanBetaWeigths'
            output_values{i} = mean(partition_betaWeights,2);
        case 'testIdx'
            output_values{i} = testIdx;
        case 'intercept'
            output_values{i} = intercept;
    end
end
end
