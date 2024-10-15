function [reduced_data, coeff] = pca_pipeline(data, ndims, varargin)

warning('off', 'all');

opts = varargin{1};
% Ensure that the data is nTrials x nDimensions (PCA expects observations in rows)
data = data';  % Transpose to get nTrials x nDimensions

% Perform PCA on the transposed data
[coeff, score, ~] = pca(data);  % score is the transformed data, coeff are the principal components

% Select the top 'newnDimensions' principal components
reduced_data = score(:, 1:ndims);

% Transpose back to get the reduced data in newnDimensions x nTrials
reduced_data = reduced_data';
warning('on', 'all');
end


