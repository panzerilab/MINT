function output_values = pca_wrapper(input, varargin)
% *function output_values = pca_wrapper(input, outputs, opts)*
%
% Usage:
%   output_values = pca_wrapper(input, varargin)
%
% Inputs:
%   - input: A cell array containing:
%       - input{1}: Feature matrix (data), where rows represent features and columns
%         represent observations.
%
%   - varargin: Optional arguments, which may include:
%       - outputs: A cell array specifying which outputs to return (default is {'all'}).
%       - opts: A structure containing options, including:
%           - standardize: Boolean indicating whether to standardize the data (default is true).
%           - numComponents: (Optional) Number of principal components to return. 
%           - explainedVariance: (Optional) Desired amount of variance to explain (default is 0.95).
%             This option is used to determine how many components to return if numComponents is NaN.
%           - NaN_handling: Specifies how NaN values should be handled in the data.                         
%           Note: you can either specify numComponents or explained Variance. Specifying both will give an error
%          
% Outputs:
%   - output_values: A cell array containing the requested output values, which may include:
%       - 'coeff': Principal component coefficients (loadings).
%       - 'score': Principal component scores.
%       - 'latent': Eigenvalues of the covariance matrix (explained variance).
%       - 'explained': Percentage of variance explained by each component.
%       - 'mean': Mean of the original data.

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
% Step 1: Check input, Check OutputList, Fill missing opts with default values  %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if nargin < 1
    error('PCA:notEnoughInput', 'Not enough input arguments. See help pca_wrapper for usage info.');
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
data = input{1}';

if ~isfield(opts, 'standardize')
    opts.standardize = true;
end

if ~isfield(opts, 'explainedVariance') && ~isfield(opts, 'numComponents')
    opts.explainedVariance = 0.95;
    opts.numComponents = NaN;
elseif ~isfield(opts, 'explainedVariance') && isfield(opts, 'numComponents')
    opts.explainedVariance = NaN;
elseif isfield(opts, 'explainedVariance') && ~isfield(opts, 'numComponents')
    opts.numComponents = NaN;
elseif ~isnan(opts.explainedVariance) && ~isnan(opts.numComponents)
    error('PCA:InvalidOptions', 'You can not set the Number of Components and the explained Variance at the same time. See help pca_wrapper for usage info.'); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                           Step 2: Perform PCA                                 %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if isnan(opts.numComponents)
    [coeff, score, latent, ~, explained, mu] = pca(data, 'Centered', opts.standardize);
else
    [coeff, score, latent, ~, explained, mu] = pca(data, 'Centered', opts.standardize, 'NumComponents', opts.numComponents);
end

% Select number of components to reach explained variance threshold
if ~isnan(opts.explainedVariance)
    cumExplained = cumsum(explained);
    numComponents = find(cumExplained >= opts.explainedVariance * 100, 1);
    coeff = coeff(:, 1:numComponents);
    score = score(:, 1:numComponents);
    latent = latent(1:numComponents);
    explained = explained(1:numComponents);
elseif ~isnan(opts.numComponents)
    numComponents = opts.numComponents;
    coeff = coeff(:, 1:numComponents);
    score = score(:, 1:numComponents);
    latent = latent(1:numComponents);
    explained = explained(1:numComponents);
end



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%                                Step 4: Outputs                                %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

possibleOutputs = {'coeff', 'score', 'latent', 'explained', 'mean'};
if ismember('all', outputs)
    outputs = possibleOutputs;
end
[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    error('PCA:invalidOutput', 'Invalid Outputs: %s', strjoin(nonMembers, ', '));
end

% Prepare the output
output_values = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'coeff'
            output_values{i} = coeff;
        case 'score'
            output_values{i} = score;
        case 'latent'
            output_values{i} = latent;
        case 'explained'
            output_values{i} = explained;
        case 'mean'
            output_values{i} = mu;
    end
end

end