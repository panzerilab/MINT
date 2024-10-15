function [predictedLabels, labelsTest, dataTest, posteriorProbs, varargout] = buildML_v4(data,labels,testSet,opts)

%%% *function [predictedLabels, labelsTest, dataTest, varargout] = buildML(data, labels, testSet, opts)*
%%%
%%% ### Description
%%% buildML trains a Machine Learning model with the chosen cross-validation algorithm on a specific set of data points and corresponding labels.
%%%
%%% ### Inputs:
%%% - *data*: n_samples x n_features predictor array.
%%% - *labels*: n_samples x 1 array of labels to be learnt.
%%% - *testSet*: fraction of trials that will be used for testing. Must be between 0 and 1. It can also be a vector with indices to select test trials from data.
%%% - *opts*: options used to build the model (see further notes).
%%%
%%% ### Outputs:
%%% - *predictedLabels*: predicted labels for samples in the held out set.
%%% - *labelsTest*: true labels for samples in the held out set.
%%% - *dataTest*: data used for the held out set.
%%% - *posteriorProbs* (only with SVM): probabilities of each label given the trial.
%%%
%%% ### Further notes:
%%% The *opts* structure can have the following fields:
%%%
%%% | field                                | description                                                | allowed values                                                                                                                                         |    default        |
%%% |--------------------------------------|------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------------|-------------------|
%%% | opts.algorithm                       | specifies the algorithm used                               | `'GLM'` (Generalized Linear Model)<br>`'linear_SVM'` (linear Support Vector Machine)<br>`'RBF_SVM'` (Radial Basis Function Support Vector Machine)     |   no default      |
%%% | opts.cv_type                         | Cross Validation type                                      | `'KFold'`<br>`'LeaveOneOut'`<br>                                                                                                                       |   no default      |
%%% | opts.K                               | number of folds (needed only if opts.cv_type = `'KFold'`)  | int > 1                                                                                                                                                |   no default      |
%%% | opts.optimize_params                 | wheter to keep default values or optimize them             | boolean                                                                                                                                                |   false           |
%%% | opts.optim_reps                      | number of parameters used for hyperparameter optim.        | int > 1                                                                                                                                                |   10              |
%%% | opts.hp_C                            | hyperparam C of the SVM                                    | float > 0                                                                                                                                              |   1               |
%%% | opts.hp_gamma                        | hyperparam gamma of the SVM                                | float > 0                                                                                                                                              |   1/num. features |
%
%  This source code is part of:
%  NIT - Neuroscience Information Toolbox
%  Copyright (C) 2023  Roberto Maffulli, Miguel Angel Casal Santiago, Marco
%  Celotto
%
%  This program is free software: you can redistribute it and/or modify
%  it under the terms of the GNU General Public License as published by
%  the Free Software Foundation, either version 3 of the License, or
%  (at your option) any later version.
%
%  This program is distributed in the hope that it will be useful,
%  but WITHOUT ANY WARRANTY; without even the implied warranty of
%  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%  GNU General Public License for more details.
% 
%  You should have received a copy of the GNU General Public License
%  along with this program.  If not, see <https://www.gnu.org/licenses/>.


%% Parameters checks
if any(isnan(data))
    error("data contains NaNs. Aborting.")
end
if any(isnan(labels))
    error("labels contains NaNs. Aborting.")
end
if any(isnan(testSet))
    error("testSet contains NaNs. Aborting.")
end
if ~isfield(opts,'verbose')
    opts.verbose = 1;
end
if ~isfield(opts,'optimize_params')
    opts.optimize_params = false;
end

assert(strcmp(opts.algorithm,'linear_SVM') ||...
    strcmp(opts.algorithm,'RBF_SVM') ||...
    strcmp(opts.algorithm, 'GLM'),...
    ['''algorithm'' argument can be only ''linear_SVM'', ''RBF_SVM'' or ''GLM''. ',...
    'Specified value ''%s'''], opts.algorithm);

GLM_flag = strcmp(opts.algorithm,'GLM');

if any(strcmp(opts.algorithm,["linear_SVM","RBF_SVM"]))
    SVM_flag = true;
    if strcmp(opts.algorithm,"linear_SVM")
        svm_family = 0;
    else
        svm_family = 2;
    end
end

if opts.optimize_params % options controls if we want to do hyperparams optimization
    
    % info parameters used for hyperparameters selection (optimal
    % hyperparameter is the one maximizing information between real and
    % predicted labels on the training set)
    info_opts.method = 'dr'; 
    info_opts.bias = 'naive'; 
    info_opts.bin_methodX = 'none'; 
    info_opts.bin_methodY = 'none'; 
    info_opts.verbose = 0; 
    
    assert(strcmp(opts.cv_type,'KFold') || strcmp(opts.cv_type,'LeaveOneOut'),...
        'The chosen crossvalidation method is not acceptable')
    
    if ~isfield(opts, 'optim_reps')
        opts.optim_reps = 10;
    end
    
    if  strcmp(opts.cv_type,'KFold')
        assert(isfield(opts,'K'),...                     
            ['With KFold crossvalidation you need to provide',...
            ' the amount of subsets (k)'])
    end
    
    if opts.verbose 
        if  isfield(opts,'hp_C') || isfield(opts,'hp_gamma')
            disp('opts.hp_C and opts.hp_gamma will be ignored')
        end
    

        if isfield(opts,'K') && strcmp(opts.cv_type,'LeaveOneOut')
            warning('''opts.K'' argument will be ignored')
        end
    end
else 
    if ~isfield(opts, 'hp_C') % Setting default hyperparam C to 1
        opts.hp_C = 1;
    end
    varargout{1} = opts.hp_C;
    if svm_family == 2
        if ~isfield(opts, 'hp_gamma') % Setting default hyperparam gamma to 1/n_features
            opts.hp_gamma = 1/size(data,2);
            varargout{1} = [opts.hp_C, opts.hp_gamma];
        end
    end
    if opts.verbose 
        if isfield(opts,'cv_type')
            warning('No hyperparameters optimization selected, opts.cv_type will be ignored')
        end
    end
end
    
%% Select training and test set
nTrials = length(labels);
nTrialsTest = length(testSet);

nTrialsTrain = nTrials-nTrialsTest;

if nTrialsTrain==0 || nTrialsTest==0                                    %check if data is correctly splitted for train and test
    error(['The number of trials for test or train cannot be zero.'...
        ' Modify your test fraction'])
end

train_test_partition.test = testSet;
train_test_partition.training = setdiff(1:size(data,1),testSet);

% Training and test set for main computation
dataTrain = data(train_test_partition.training',:);                     %labels for train
labelsTrain = labels(train_test_partition.training');                   %trials for train
dataTest = data(train_test_partition.test',:);                          %labels for test
labelsTest = labels(train_test_partition.test');                        %trials for test

% Training and test sets for hyperparameters optimization
if opts.optimize_params % cross validation folds on the training set, to do hyperparam optimization
    if strcmp(opts.cv_type,'KFold')
        cvPartition = cvpartition(labelsTrain,'KFold',opts.K);              %split data for kFold crossvalidation
    elseif strcmp(opts.cv_type,'LeaveOneOut')
        cvPartition = cvpartition(labelsTrain,'LeaveOut');                  %split data for LeaveOneOut crossvalidation
    end

end

%% Run SVM 
if SVM_flag                                                             % if SVM is chosen  
    if opts.optimize_params                                             % if we want to do hyperparameter optimization on the training set
        [repInfo,svmvarout] = optimize_hyperparams(dataTrain, labelsTrain, cvPartition, svm_family, opts, info_opts);

        % select optimal hyperparam for linear SVM
        if svm_family == 0
            [~,c_idx] = max(repInfo);
            c_opt = svmvarout{1}(c_idx); % svmvarout{1} = C values spanned during hyperparam optimization
            varargout{1} = c_opt;
            opts.hp_C = c_opt;
            
        % select optimal hyperparam for RBF kernel SVM
        elseif svm_family == 2
             [~,param_idx] = max(repInfo(:));
             [c_idx,g_idx] = ind2sub(size(repInfo),param_idx); 
             c_opt = svmvarout{1}(c_idx); % svmvarout{1} = C values spanned during hyperparam optimization
             gamma_opt = svmvarout{2}(g_idx); % svmvarout{2} = gamma values spanned during hyperparam optimization
             varargout{1} = [c_opt, gamma_opt];
             opts.hp_C = c_opt;
             opts.hp_gamma = gamma_opt;
        end
        
    else                                                                    % we want to use the default SVM hyperparameters
        % No hyperparameter optimization
    end
    % Train final SVM model on this fold 
    if svm_family == 0
        algorithm_opt = svmtrain(labelsTrain, dataTrain,...
                 sprintf('-q -c %f -t %d -b 1', opts.hp_C,...
                 svm_family));
    elseif svm_family == 2
        algorithm_opt = svmtrain(labelsTrain, dataTrain,...
                 sprintf('-q -c %f -g %f -t %d -b 1', opts.hp_C, opts.hp_gamma,...
                 svm_family));
    end

   % Test final SVM model on this fold     
    [predictedLabels, ~, posteriorProbs] = ...
             svmpredict(labelsTest, dataTest, algorithm_opt, '-q -b 1'); 
   
        
elseif GLM_flag                                                         %if GLM chosen 
    foldId = zeros(size(labelsTrain));
    
    for i = 1:setsCV                                                    %for each set of CV
        foldId = foldId + i*cvPartition.test(i);
    end
    
    CVfit = cvglmnet(dataTrain,labelsTrain,...
        'multinomial',[],'class',[],foldId);
    predictedLabels = cvglmnetPredict(CVfit,dataTest,...
        'lambda_min','class');
end
end


function [repInfo,svmvarout] = optimize_hyperparams(dataTrain, labelsTrain, cvPartition, svm_family, opts, info_opts)

setsCV = cvPartition.NumTestSets;                                       %crossvalidation sets

c_vals=logspace(-2,2,opts.optim_reps);                          % vector of C parameters to try (logaritmically equispaced between 10^-2 and 100)
svmvarout{1} = c_vals;

repInfo = zeros(1,opts.optim_reps);                             % population information across repetitions
if svm_family == 2 % RBF kernel, also gamma hyperparam
    g_vals=logspace(-5,1,opts.optim_reps);                      % vector of gamma parameters to try (logaritmically equispaced between 10^-5 and 10)
    svmvarout{2}=g_vals;
    
    repInfo = zeros(opts.optim_reps,opts.optim_reps);
end
for repIdx1 = 1:opts.optim_reps % loop over number of optimization parameters
    optimPredictedLabels = nan(1,numel(labelsTrain));
    if svm_family == 0 % linear, just one hyperparam (C)
        for i = 1:setsCV                                                    %for each set of CV
            testIdx = find(cvPartition.test(i));                            %trials for test CV
            trainIdx = find(cvPartition.training(i));
            
            algorithm_opt = svmtrain(labelsTrain(trainIdx), dataTrain(trainIdx,:),...
                sprintf('-q -c %f -t %d -b 1', c_vals(repIdx1), ...
                svm_family)); 
            [optimPredictedLabels(testIdx), ~, posteriorProbs] = ...
                svmpredict(labelsTrain(testIdx), dataTrain(testIdx,:), algorithm_opt, '-q -b 1');
        end
        repInfo(repIdx1) = cell2mat(information(labelsTrain',optimPredictedLabels,info_opts,{'I'}));
    
    elseif svm_family == 2 % RBF kernel, 2 hyperparameters (C and gamma)
        for repIdx2 = 1:opts.optim_reps % loop over gamma
            optimPredictedLabels = nan(1,numel(labelsTrain));
            for i = 1:setsCV                                                    %for each set of CV
                testIdx = find(cvPartition.test(i));                            %trials for test CV
                trainIdx = find(cvPartition.training(i));
                
                algorithm_opt = svmtrain(labelsTrain(trainIdx), dataTrain(trainIdx,:),...
                    sprintf('-q -c %f -g %f -t %d -b 1', c_vals(repIdx1), g_vals(repIdx2), ...
                    svm_family)); 
                [optimPredictedLabels(testIdx), ~, posteriorProbs] = ...
                    svmpredict(labelsTrain(testIdx), dataTrain(testIdx,:), algorithm_opt, '-q -b 1');
            end
            repInfo(repIdx1,repIdx2) = cell2mat(information(labelsTrain',optimPredictedLabels,info_opts,{'I'}));
        end
    end
end

end
