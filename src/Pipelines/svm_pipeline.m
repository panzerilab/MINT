function [predictedLabels, labelsTest, dataTest, posteriorProbs, optimHyp, betaWeights, intercept] = svm_pipeline(data, labels, testSet, opts)
%%% svm_pipeline: Train and evaluate Support Vector Machine (SVM) models with linear or RBF kernel using LIBSVM or FITCSVM.
%%%   This function allows training and testing SVM models with either linear or RBF kernel, using either LIBSVM or FITCSVM as the underlying SVM implementation.
%%%   It includes options for hyperparameter optimization. The function outputs predicted labels,
%%%   true labels, test data, and posterior probabilities and, if hyperparameter optimization is enabled, optimal
%%%   hyperparameters as additional output.
%%%
%%% Input:
%%%   - data: Training data matrix (each row represents a sample, each column a feature)
%%%   - labels: True labels corresponding to the training data samples
%%%   - testSet: Indices of samples to be used as the test set
%%%   - opts: Structure with optional parameters for configuring the SVM model and optimization process
%%%       - optimize_params: Flag for hyperparameter optimization (default is false)
%%%       - cv_type: Cross-validation type, 'KFold' or 'LeaveOneOut'
%%%       - K: Number of folds for KFold cross-validation (required if cv_type is 'KFold')
%%%       - optim_reps: Number of repetitions for hyperparameter optimization (default is 10)
%%%       - svm_family: Type of SVM kernel, 'linear' or 'RBF' for Radial Basis Function
%%%       - hp_C: Hyperparameter C for SVM regularization (default is 1)
%%%       - hp_gamma: Hyperparameter gamma for RBF kernel (default is 1/n_features)
%%%       - libsvm: Flag for using LIBSVM (true) or FITCSVM (false) (default is false)
%%%
%%% Output:
%%%   - predictedLabels: Predicted labels for the test set
%%%   - labelsTest: True labels for the test set
%%%   - dataTest: Test data matrix corresponding to the test set
%%%   - posteriorProbs: Posterior probabilities of class membership for the test set
%%%   - optimHyp: If hyperparameter optimization is performed, optimHyp contains optimal hyperparameters


warning('off', 'all');

%% Check input parameters
optimHyp = struct();

% Validate input parameters and set default values if necessary
if any(isnan(data(:))) || any(isnan(labels)) || any(isnan(testSet(:)))
    error("One or more variables contain NaNs. Aborting.")
end

if ~isfield(opts,'optimize_params')
    opts.optimize_params = false;
end
if ~isfield(opts,'libsvm')
    opts.libsvm = false;
end
if opts.optimize_params
    if opts.libsvm == true
       error("Hyperparameter optimization is not supported for libsvm in this pipeline. Please set the 'libsvm' parameter to false to enable hyperparameter optimization.")
    end 

    if ~isfield(opts,'gamma_range')
        gamma_range = [1e-3, 1e3];
    else
        gamma_range = opts.gamma_range;
    end

    if ~isfield(opts,'C_range')
        C_range =[ 1e-3, 1e3];
    else
        C_range = opts.C_range;
    end

    if ~isfield(opts,'parallel_optim')
        opts.parallel_optim = false;
    end
else
    if ~isfield(opts, 'hp_C')
        opts.hp_C = 1;
    end
    varargout{1} = opts.hp_C;
    if ~isfield(opts, 'Standardize')
        opts.Standardize = true;
    end 
    if  strcmp(opts.svm_family, 'RBF')
        if ~isfield(opts, 'hp_gamma')
            opts.hp_gamma = 1/size(data,2);
            optimHyp.C = opts.hp_C;
            optimHyp.gamma = opts.hp_gamma;
        end
    end
end

if isfield(opts,'K') && strcmp(opts.cv_type,'LeaveOneOut')
    warning('''opts.K'' argument will be ignored')
end

if  strcmp(opts.cv_type,'KFold')
    assert(isfield(opts,'K'),...
        'With KFold crossvalidation you need to provide the amount of subsets (k)')
end

if ~(strcmp(opts.cv_type, 'KFold') || strcmp(opts.cv_type, 'LeaveOneOut'))
    error('The chosen cross-validation method is not acceptable');
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


%% Check if labels are binary

unique_values = unique(labels);

if size(unique_values, 1) > 2
    opts.binary = false;
else
    opts.binary = true;
end

%% Linear SVM

% Train and evaluate Linear SVM
if strcmp(opts.svm_family, 'linear')
    if opts.binary == true
        if opts.libsvm == true 
            svmModel = svmtrain(labelsTrain, dataTrain,sprintf('-q -c %f -t %d -b 1', opts.hp_C,0));
            betaWeights = svmModel.SVs' * svmModel.sv_coef;
            intercept = -svmModel.rho;
            [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
        else
            if opts.optimize_params
                hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
                                       optimizableVariable('Standardize', {'true', 'false'})];
                opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus', 'UseParallel', opts.parallel_optim,  'MaxObjectiveEvaluations', 30, 'ShowPlots', false, 'Verbose', 0, opts.cv_type, opts.K);
                svmModel = fitcsvm(dataTrain, labelsTrain, 'KernelFunction', 'linear', 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt);
                optimHyp.BoxConstraint = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.BoxConstraint; 
                optimHyp.Standardize = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.Standardize;
                svmModel = fitPosterior(svmModel);
                betaWeights = svmModel.Beta;
                intercept = svmModel.Bias;
                [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
            else
                svmModel = fitcsvm(dataTrain, labelsTrain, 'KernelFunction', 'linear', 'BoxConstraint', opts.hp_C, 'Standardize', true);
                svmModel = fitPosterior(svmModel);
                betaWeights = svmModel.Beta;
                intercept = svmModel.Bias;
                [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
            end
        end
    else
        if opts.libsvm == true
            svmModel = svmtrain(labelsTrain, dataTrain, sprintf('-q -s 0 -c %f -t %d -b 1', opts.hp_C, 0));
            [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
        else
            if opts.optimize_params
                hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
                                       optimizableVariable('Standardize', {'true', 'false'})];
                opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus', 'UseParallel', opts.parallel_optim, 'MaxObjectiveEvaluations', 30, 'ShowPlots', false, 'Verbose', 0, opts.cv_type, opts.K);
                svmTemplate = templateSVM('KernelFunction','linear');
                svmModel = fitcecoc(dataTrain, labelsTrain, 'Learners', svmTemplate, 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt, 'Coding', 'onevsone');
                optimHyp.BoxConstraint = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.BoxConstraint; 
                optimHyp.Standardize = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.Standardize;
                binaryLearners = svmModel.BinaryLearners;
                betaWeights = [];
                intercept = [];
                for i = 1:numel(binaryLearners)
                    binaryModel = binaryLearners{i};
                    if ~isempty(binaryModel)
                        betaWeights = [betaWeights, binaryModel.Beta]; 
                        intercept = [intercept, binaryModel.Bias];
                    end
                end
                [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
            else
                svmTemplate = templateSVM('KernelFunction','linear', 'BoxConstraint', opts.hp_C, 'Standardize', true);
                svmModel = fitcecoc(dataTrain, labelsTrain, 'Learners',  svmTemplate, 'Coding', 'onevsone');
                binaryLearners = svmModel.BinaryLearners;
                betaWeights = [];
                intercept = [];
                for i = 1:numel(binaryLearners)
                    binaryModel = binaryLearners{i};
                    if ~isempty(binaryModel)
                        betaWeights = [betaWeights, binaryModel.Beta]; 
                        intercept = [intercept, binaryModel.Bias];
                    end
                end
                [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
            end
        end
    end
end



%% RBF SVM 
% Train and evaluate RBF SVM
if strcmp(opts.svm_family, 'RBF')
    if opts.binary == true
        if opts.libsvm == true
            svmModel = svmtrain(labelsTrain, dataTrain,sprintf('-q -c %f -g %f -t %d -b 1', opts.hp_C, opts.hp_gamma,2));
            betaWeights = svmModel.SVs' * svmModel.sv_coef;
            intercept = -svmModel.rho;
            [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
        else
            if opts.optimize_params
                hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
                                       optimizableVariable('KernelScale', gamma_range, 'Transform', 'log'), ...
                                       optimizableVariable('Standardize', {'true', 'false'})];
                opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus', 'UseParallel', opts.parallel_optim, 'MaxObjectiveEvaluations', 30, 'ShowPlots', false, 'Verbose', 0, opts.cv_type, opts.K);
                svmModel = fitcsvm(dataTrain, labelsTrain, 'KernelFunction', 'rbf', 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt);
                optimHyp.BoxConstraint = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.BoxConstraint; 
                optimHyp.Standardize = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.Standardize;
                optimHyp.KernelScale = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.KernelScale;                
                svmModel = fitPosterior(svmModel);
                betaWeights = svmModel.Beta;
                intercept = svmModel.Bias;
                [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
            else
                svmModel = fitcsvm(dataTrain, labelsTrain, 'KernelFunction', 'rbf', 'KernelScale', (1/(2*opts.hp_gamma)),'BoxConstraint', opts.hp_C, 'Standardize', opts.Standardize);
                svmModel = fitPosterior(svmModel);
                betaWeights = svmModel.Beta;
                intercept = svmModel.Bias;
                [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
            end
        end
    else
        if opts.libsvm == true
            svmModel = svmtrain(labelsTrain, dataTrain,sprintf('-q -s 0 -c %f -g %f -t %d -b 1', opts.hp_C, opts.hp_gamma,2));
            [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
        else
            if opts.optimize_params
                hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
                                       optimizableVariable('KernelScale', gamma_range, 'Transform', 'log'),...
                                       optimizableVariable('Standardize', {'true', 'false'})];
                opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus','UseParallel', opts.parallel_optim, 'MaxObjectiveEvaluations', 30, 'ShowPlots', false, 'Verbose', 0, opts.cv_type, opts.K);
                svmTemplate = templateSVM('KernelFunction','rbf');
                svmModel = fitcecoc(dataTrain, labelsTrain, 'Learners', svmTemplate, 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt, 'Coding', 'onevsone');
                optimHyp.BoxConstraint = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.BoxConstraint; 
                optimHyp.Standardize = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.Standardize;
                optimHyp.KernelScale = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.KernelScale;
                binaryLearners = svmModel.BinaryLearners;
                betaWeights = [];
                intercept = [];
                for i = 1:numel(binaryLearners)
                    binaryModel = binaryLearners{i};
                    if ~isempty(binaryModel)
                        betaWeights = [betaWeights, binaryModel.Beta]; 
                        intercept = [intercept, binaryModel.Bias];
                    end
                end
                [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
            else
                svmTemplate = templateSVM('KernelFunction','rbf','KernelScale', (1/(2*opts.hp_gamma)), 'BoxConstraint', opts.hp_C, 'Standardize', opts.Standardize);
                svmModel = fitcecoc(dataTrain, labelsTrain, 'Learners',  svmTemplate, 'Coding', 'onevsone');
                binaryLearners = svmModel.BinaryLearners;
                betaWeights = [];
                intercept = [];
                for i = 1:numel(binaryLearners)
                    binaryModel = binaryLearners{i};
                    if ~isempty(binaryModel)
                        betaWeights = [betaWeights, binaryModel.Beta]; 
                        intercept = [intercept, binaryModel.Bias];
                    end
                end
                [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
            end
        end
    end
end
warning('on', 'all');
end


