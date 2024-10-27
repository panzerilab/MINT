function output_values = svm_wrapper(input,varargin)
% *function output_values = svm_wrapper(input, outputs, opts)*
% SVM - Support Vector Machine Implementation
%
% This function implements a Support Vector Machine (SVM) for binary and multi-class classification
% using various kernel functions. It includes options for hyperparameter optimization and cross-validation.
%
% Usage:
%   output_values = SVM(input, varargin)
%
% Inputs:
%   - input: A cell array where:
%       - input{1}: Feature matrix (data), where rows represent observations and columns represent features.
%       - input{2}: Label vector (labels) corresponding to the observations in the feature matrix.
%
%   - varargin: Optional arguments, which can include:
%       - outputs: A cell array specifying which outputs to return (default is {'all'}). 
%       - opts: A structure containing options that include:
%           - libsvm: Boolean indicating whether to use LIBSVM (default false).
%           - optimize_params: A structure with optimization options for hyperparameters (default is empty).
%           - svm_family: Type of SVM to use, options include 'linear' and 'RBF' (default 'linear').
%           - cv_type: Type of cross-validation (default 'KFold').
%           - K_folds: Number of folds for K-fold cross-validation (default 5).
%           - hp_C: Regularization parameter for SVM (default 1).
%           - hp_gamma: Kernel coefficient for 'RBF' SVM (default 1/number of observations for RBF).
%           - Standardize: Boolean indicating whether to standardize features (default true).
%           - optim_opts: A structure with options for hyperparameter optimization which can include:
%               - optimize: Boolean indicating whether to optimize hyperparameters (default false).
%               - gamma_range: Range for gamma values during optimization (default [1e-3, 1e3]).
%               - C_range: Range for C values during optimization (default [1e-3, 1e3]).
%               - parallel: Boolean indicating whether to run optimization in parallel (default false).
%               - MaxIter: Maximum iterations for optimization (default 50).
%               - optim_cv: Cross-validation method for optimization (default {'Resubstitution'}).
%
% Outputs:
%   - output_values: A cell array containing the requested output values, which may include:
%       - 'labels': Predicted class labels for the test data.
%       - 'posteriorProbs': Posterior probabilities for each class.
%       - 'optimized_hyperparameters': Hyperparameters used for the SVM.
%       - 'betaWeights': Weights of the SVM model (only for linear SVM).
%       - 'intercept': Intercept term of the SVM model.
%       - 'mean_betaWeights': Mean of the beta weights across all folds (only for linear SVM).
%       - 'testIdx': Indices of the test set used in each fold to fit the SVM.
%
% Note: Hyperparameter optimization is supported for both linear and RBF SVMs, but certain configurations
% (e.g., LIBSVM) may restrict options for optimization. Ensure appropriate parameters are set before calling
% the function to avoid errors.
%
% Example:
%   input = {dataMatrix, labelVector};
%   opts = struct('svm_family', 'RBF', 'libsvm', false);
%   outputs = {'labels', 'posteriorProbs'};
%   output_values = SVM(input, outputs, opts);

% Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
% This file is part of MINT. 
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step 1: Check input, Check OutputList, Fill missing opts with default values %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
defaultOpts.libsvm = false;
defaultOpts.optimize_params = struct();
defaultOpts.svm_family = 'linear';
defaultOpts.cv_type = 'KFold';
defaultOpts.K_folds = 5;
defaultOpts.NaN_handling = 'error';

if nargin < 1
    msg = 'not enough input arguments. See `help SVM` for usage info';
    error('SVM:notEnoughInput', msg);
elseif nargin == 1
    opts = defaultOpts;
    outputs = {'all'};
elseif nargin == 2
    if iscell(varargin{1})
        outputs = varargin{1};
        opts = defaultOpts;
    elseif isstruct(varargin{1})
        opts = varargin{1};
        outputs = {'all'};
    end
elseif nargin == 3
    opts = varargin{2};
    outputs = varargin{1};
end
if ~isfield(opts, 'NaN_handling')
    opts.NaN_handling = 'error';
end
input = nan_method(input, opts.NaN_handling);
% Load the data
data = input{1};
labels = input{2};

% check opts
if ~isfield(opts,'svm_family')
    opts.svm_family = 'linear';
end
if ~isfield(opts,'optim_opts')
    optimization_opts.optimize = false;
else
    optimization_opts = opts.optim_opts;
    if ~isfield( optimization_opts, 'optimze')
        optimization_opts.optimize = true;
    end
end

if ~isfield(opts,'libsvm')
    opts.libsvm = false;
end

if optimization_opts.optimize == true
    hp_C = NaN;
    hp_gamma = NaN;
    opts.Standardize = NaN;
    if opts.libsvm == true
        error("Hyperparameter optimization is not supported for libsvm in this pipeline. Please set the 'libsvm' parameter to false to enable hyperparameter optimization.")
    end
    if ~isfield(optimization_opts,'gamma_range')
        optimization_opts.gamma_range = [1e-3, 1e3];
    end

    if ~isfield(optimization_opts,'C_range')
        optimization_opts.C_range =[ 1e-3, 1e3];
    end
    if ~isfield(optimization_opts,'parallel')
        optimization_opts.parallel = false;
    end
    if ~isfield(optimization_opts,'MaxIter')
        optimization_opts.MaxIter = 50;
    end

    if ~isfield(optimization_opts,'optim_cv')
        optimization_opts.optim_cv = {'Resubstitution'};
        warning('No cross-validation method specified for the hyperparameter optimization. Using default method: ''Resubstitution''. The SVM will use the entire dataset to optimize the hyperparameters.');
    end
    optim_cv_method = string(optimization_opts.optim_cv{1});
    valid_methods = ["KFold", "HoldOut", "LeaveMOut", "LeaveOut", "Resubstitution"];
    if ~ismember(optim_cv_method, valid_methods)
        error('SVM:The chosen cross-validation method for hyperparameter optimization is not acceptable');
    end
    if strcmp(optim_cv_method,'KFold')
        if length(optimization_opts.optim_cv) > 1
            optim_cv_partition =  optimization_opts.optim_cv{2};
        else
            optim_cv_partition =  5;
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
else
    if ~isfield(opts, 'hp_C')
        hp_C = 1;
    else
        hp_C = opts.hp_C;
    end
    if ~isfield(opts, 'Standardize')
        opts.Standardize = true;
    end
    if  strcmp(opts.svm_family, 'RBF')
        if ~isfield(opts, 'hp_gamma')
            hp_gamma = 1/size(data,1);
        else
            hp_gamma = opts.hp_gamma;
        end
    end
    optim_cvPartition = NaN;
end

if ~isfield(opts,'cv')
    opts.cv = {'Resubstitution'};
    warning('No cross-validation method specified. Using default method: ''Resubstitution''. The SVM will be fitted on the entire dataset without cross-validation.');
end

if iscell(opts.cv)
    cv_method =  opts.cv{1};
    valid_methods = ["KFold", "HoldOut", "LeaveMOut", "LeaveOut", "Resubstitution"];
    if ~ismember(cv_method, valid_methods)
        error('SVM:The chosen cross-validation method for hyperparameter optimization is not acceptable');
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
    error('Invalid input for cross-validation. The cross-validation option must be provided as a cell array. For ''KFold'', specify the number of folds as {''KFold'', numFolds}, e.g., {''KFold'', 5}. The other valid option is {''LeaveOneOut''}.');
end

% Check outputs
possibleOutputs = {'labels', 'posteriorProbs', 'optimized_hyperparameters', 'betaWeights', 'intercept', 'mean_betaWeigths', 'testIdx'};
if ismember('all', outputs)
    if strcmp(opts.svm_family, 'RBF')
        outputs = {'labels', 'posteriorProbs', 'optimized_hyperparameters', 'intercept'};
    else         
        outputs = {'labels', 'posteriorProbs', 'optimized_hyperparameters', 'betaWeights', 'intercept', 'mean_betaWeigths', 'testIdx'};
    end 
end 
if strcmp(opts.svm_family, 'RBF')
    invalidFields = {'betaWeights', 'mean_betaWeigths'};
    [isInvalid, invalidIdx] = ismember(invalidFields, outputs);
    if any(isInvalid)
        outputs(invalidIdx(isInvalid)) = [];
        warning('SVM:invalidOutput', 'betaWeights and mean_betaWeights are not supported for RBF SVM and have been removed from the outputs.');
    end
end
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    msg = sprintf('Invalid Outputs: %s', strjoin(nonMembers, ', '));
    error('SVM:invalidOutput', msg);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Step 2: Binary or Non Binary?                       %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
unique_values = unique(labels);
if size(unique_values, 2) > 2
    isbinary = false;
else
    isbinary = true;  
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%              Step 3: Select training and test set and fit svm                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
predictedLabels = zeros(length(labels), 1);
posteriorProbability = zeros(length(labels),length(unique_values));
if isbinary
    intercept = zeros(1, cvPartition.NumTestSets);
    if strcmp(opts.svm_family, 'linear')
        partition_betaWeights = zeros(size(input{1}, 1),cvPartition.NumTestSets);
    end
else
    intercept = struct();
    partition_betaWeights = struct();
end
optimHyp = struct();
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
    dataTrain = data(:,train_test_partition.training);
    labelsTrain = labels(train_test_partition.training);
    dataTest = data(:,train_test_partition.test);
    labelsTest = labels(train_test_partition.test);
    testIdx{partition} = train_test_partition.test;

    if optimization_opts.optimize == true
        if isnan(cv_partition)
            optim_cvPartition =  cvpartition(labelsTrain, optim_cv_method);
        else 
            optim_cvPartition =  cvpartition(labelsTrain, optim_cv_method, cv_partition);
        end 
    else 
        optim_cvPartition = NaN;
    end 
    % fit SVM for partition
    if strcmp(opts.svm_family, 'linear') % linear SVM
        if opts.libsvm
            [predictedLabels(test_idx), posteriorProbability(test_idx,:), partition_betaWeights.partition{partition}, intercept.partition{partition}] = fitLinearSVM_libSVM(labelsTrain', labelsTest', dataTrain', dataTest', hp_C, isbinary);
        else
            [predictedLabels(test_idx), posteriorProbability(test_idx,:), partition_betaWeights_tmp, intercept_tmp, optimHyp.partition{partition}] = fitLinearSVM(labelsTrain', dataTrain', dataTest', hp_C, optimization_opts, optim_cvPartition, opts.Standardize, isbinary);
        end
        if isbinary
            partition_betaWeights(:,partition) =  partition_betaWeights_tmp;
        else 
            partition_betaWeights.partition{partition} =  partition_betaWeights_tmp;
        end 
    elseif  strcmp(opts.svm_family, 'RBF') % RBF SVM
        if opts.libsvm
            [predictedLabels(test_idx), posteriorProbability(test_idx,:),  intercept(partition)] = fitRBFSVM_libSVM(labelsTrain', labelsTest', dataTrain', dataTest', hp_C, hp_gamma, isbinary);
        else
            [predictedLabels(test_idx), posteriorProbability(test_idx,:),  intercept.partition{partition}, optimHyp.partition{partition}] = fitRBFSVM(labelsTrain', dataTrain', dataTest', hp_C, hp_gamma, optimization_opts, optim_cvPartition, opts.Standardize, isbinary);
        end
    end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                   Step 4: Outputs                             %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
output_values = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'labels'
            output_values{i} = predictedLabels';
        case 'posteriorProbs'
            output_values{i} = posteriorProbability';
        case 'optimized_hyperparameters'
            output_values{i} = optimHyp;          
        case 'betaWeights'
            output_values{i} = partition_betaWeights;  
        case 'intercept'
            output_values{i} = intercept;
        case 'mean_betaWeigths'
            if isstruct(partition_betaWeights)
                for part = 1:cvPartition.NumTestSets
                    partition_betaWeights_all(:,:,part) = partition_betaWeights.part;                    
                end 
                output_values{i} = mean(partition_betaWeights_all,3);
            else 
                output_values{i} = mean(partition_betaWeights,2);
            end 
        case 'testIdx'
            output_values{i} = testIdx;
    end
end
end 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                               Helper Functions                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function  [predictedLabels, posteriorProbs, betaWeights, intercept] = fitLinearSVM_libSVM(labelsTrain, labelsTest, dataTrain, dataTest, hp_C, isbinary)
if isbinary
    svmModel = svmtrain(labelsTrain, dataTrain,sprintf('-q -c %f -t %d -b 1', hp_C,0));
    betaWeights = svmModel.SVs' * svmModel.sv_coef;
    intercept = -svmModel.rho;
    [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
else
    svmModel = svmtrain(labelsTrain, dataTrain, sprintf('-q -s 0 -c %f -t %d -b 1', hp_C, 0));
    betaWeights = svmModel.SVs' * svmModel.sv_coef;
    intercept = -svmModel.rho;
    [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
end
end

function  [predictedLabels, posteriorProbs, betaWeights, intercept] = fitRBFSVM_libSVM(labelsTrain, labelsTest, dataTrain, dataTest, hp_C, hp_gamma, isbinary)
if isbinary
    svmModel = svmtrain(labelsTrain, dataTrain, sprintf('-q -c %f -g %f -t %d -b 1', hp_C, hp_gamma,2));
    betaWeights = svmModel.SVs' * svmModel.sv_coef;
    intercept = -svmModel.rho;
    [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
else
    svmModel = svmtrain(labelsTrain, dataTrain,sprintf('-q -s 0 -c %f -g %f -t %d -b 1', hp_C, hp_gamma,2));
    betaWeights = svmModel.SVs' * svmModel.sv_coef;
    intercept = -svmModel.rho;
    [predictedLabels, ~, posteriorProbs] = svmpredict(labelsTest, dataTest, svmModel, '-q -b 1');
end
end



function  [predictedLabels, posteriorProbs, betaWeights, intercept, optimHyp] = fitLinearSVM(labelsTrain, dataTrain, dataTest, hp_C, optimization_opts, cvPartition, Standardize, isbinary)
optimHyp = struct();
if isbinary
    if optimization_opts.optimize
        C_range = optimization_opts.C_range;
        parallel_optim = optimization_opts.parallel;
        iterations = optimization_opts.MaxIter;
        hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
            optimizableVariable('Standardize', {'true', 'false'})];
        opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus', 'UseParallel', parallel_optim,  'MaxObjectiveEvaluations', iterations, 'ShowPlots', false, 'Verbose', 0, 'CVPartition', cvPartition);
        svmModel = fitcsvm( dataTrain, labelsTrain, 'KernelFunction', 'linear', 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt);
        optimHyp.BoxConstraint = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.BoxConstraint;
        optimHyp.Standardize = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.Standardize;
        svmModel = fitPosterior(svmModel);
        betaWeights = svmModel.Beta;
        intercept = svmModel.Bias;
        [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
    else
        svmModel = fitcsvm( dataTrain, labelsTrain, 'KernelFunction', 'linear', 'BoxConstraint', hp_C, 'Standardize', Standardize);
        svmModel = fitPosterior(svmModel);
        betaWeights = svmModel.Beta;
        intercept = svmModel.Bias;
        [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
    end
else
    if optimization_opts.optimize
        C_range = optimization_opts.C_range;
        parallel_optim = optimization_opts.parallel;
        iterations = optimization_opts.MaxIter;
        hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
            optimizableVariable('Standardize', {'true', 'false'})];
        opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus', 'UseParallel', parallel_optim, 'MaxObjectiveEvaluations', iterations, 'ShowPlots', false, 'Verbose', 0, 'CVPartition', cvPartition);
        svmTemplate = templateSVM('KernelFunction','linear');
        svmModel = fitcecoc( dataTrain, labelsTrain, 'Learners', svmTemplate, 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt, 'Coding', 'onevsone');
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
        [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
    else
        svmTemplate = templateSVM('KernelFunction','linear', 'BoxConstraint', hp_C, 'Standardize', Standardize);
        svmModel = fitcecoc( dataTrain, labelsTrain, 'Learners',  svmTemplate, 'Coding', 'onevsone');
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
        [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
    end
end
end

function  [predictedLabels, posteriorProbs, intercept, optimHyp] = fitRBFSVM(labelsTrain, dataTrain,  dataTest, hp_C, hp_gamma, optimization_opts, cvPartition, Standardize, isbinary)
optimHyp = struct();
if isbinary
    if optimization_opts.optimize
        C_range = optimization_opts.C_range;
        parallel_optim = optimization_opts.parallel;
        gamma_range = optimization_opts.gamma_range;
        iterations = optimization_opts.MaxIter;
        hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
            optimizableVariable('KernelScale', gamma_range, 'Transform', 'log'), ...
            optimizableVariable('Standardize', {'true', 'false'})];
        opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus', 'UseParallel', parallel_optim, 'MaxObjectiveEvaluations',iterations, 'ShowPlots', false, 'Verbose', 0, 'CVPartition', cvPartition);
        svmModel = fitcsvm( dataTrain, labelsTrain, 'KernelFunction', 'rbf', 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt);
        optimHyp.BoxConstraint = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.BoxConstraint;
        optimHyp.Standardize = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.Standardize;
        optimHyp.KernelScale = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.KernelScale;
        svmModel = fitPosterior(svmModel);
        intercept = svmModel.Bias;
        [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
    else
        svmModel = fitcsvm( dataTrain, labelsTrain, 'KernelFunction', 'rbf', 'KernelScale', (1/(2*hp_gamma)),'BoxConstraint', hp_C, 'Standardize', Standardize);
        svmModel = fitPosterior(svmModel);
        intercept = svmModel.Bias;
        [predictedLabels, posteriorProbs, ~] = predict(svmModel, dataTest);
    end
else
    if optimization_opts.optimize
        C_range = optimization_opts.C_range;
        parallel_optim = optimization_opts.parallel;
        gamma_range = optimization_opts.gamma_range;
        iterations = optimization_opts.MaxIter;
        hyperparameterRange = [optimizableVariable('BoxConstraint', C_range, 'Transform', 'log'),...
            optimizableVariable('KernelScale', gamma_range, 'Transform', 'log'),...
            optimizableVariable('Standardize', {'true', 'false'})];
        opt = struct('Optimizer', 'bayesopt', 'AcquisitionFunctionName', 'expected-improvement-plus','UseParallel', parallel_optim, 'MaxObjectiveEvaluations', iterations, 'ShowPlots', false, 'Verbose', 0, 'CVPartition', cvPartition);
        svmTemplate = templateSVM('KernelFunction','rbf');
        svmModel = fitcecoc( dataTrain, labelsTrain, 'Learners', svmTemplate, 'OptimizeHyperparameters', hyperparameterRange, 'HyperparameterOptimizationOptions', opt, 'Coding', 'onevsone');
        optimHyp.BoxConstraint = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.BoxConstraint;
        optimHyp.Standardize = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.Standardize;
        optimHyp.KernelScale = svmModel.HyperparameterOptimizationResults.XAtMinEstimatedObjective.KernelScale;
        binaryLearners = svmModel.BinaryLearners;
        intercept = [];
        for i = 1:numel(binaryLearners)
            binaryModel = binaryLearners{i};
            if ~isempty(binaryModel)
                intercept = [intercept, binaryModel.Bias];
            end
        end
        [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
    else
        svmTemplate = templateSVM('KernelFunction','rbf','KernelScale', (1/(2*hp_gamma)), 'BoxConstraint', hp_C, 'Standardize', Standardize);
        svmModel = fitcecoc( dataTrain, labelsTrain, 'Learners',  svmTemplate, 'Coding', 'onevsone');
        binaryLearners = svmModel.BinaryLearners;
        intercept = [];
        for i = 1:numel(binaryLearners)
            binaryModel = binaryLearners{i};
            if ~isempty(binaryModel)               
                intercept = [intercept, binaryModel.Bias];
            end
        end
        [predictedLabels, posteriorProbs] = predict(svmModel, dataTest);
    end
end
end




