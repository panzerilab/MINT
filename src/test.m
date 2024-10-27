clc, clear all;
rng(42);

nSamples = 1000;   
nFeatures = 10;   
noiseLevel = 0.05;  

rng(42);  % For reproducibility
X = randn(nSamples, nFeatures);  % Random features (normal distribution)
trueBeta = randn(nFeatures, 1);  % True coefficients
intercept = 0.5;  % True intercept

% Generate linear combination of features and add intercept
linearCombination = X * trueBeta + intercept;

% Convert linear combination to probabilities using the logistic function
probabilities = 1 ./ (1 + exp(-linearCombination));  % Sigmoid function

% Generate binary labels using these probabilities
labels = binornd(1, probabilities);  % Bernoulli distribution (0 or 1)

% Add noise to the labels (optional, can be turned off)
if noiseLevel > 0
    noisyLabels = labels;
    flipIdx = rand(size(labels)) < noiseLevel;  % Flip a percentage of the labels
    noisyLabels(flipIdx) = 1 - noisyLabels(flipIdx);  % Flip 0 to 1 and 1 to 0
else
    noisyLabels = labels;
end

opts.cv = {'KFold', 4};
opts.glmnet = false;
optimization.optim_cv = {'KFold', 2};
opts.optim_opts = optimization;
opts.distribution = 'normal';
opts.regularization = 'elasticNet';

outputs = {'all'};

output_values = glm_wrapper({X', labels},  outputs, opts);

labels_pred = output_values{1};
labelsBinary = output_values{2}';
lambda = output_values{3};
confusionMatrix = output_values{4};
betaWeights = output_values{5};
meanBetaWeigths = output_values{6};
testIdx = output_values{7};