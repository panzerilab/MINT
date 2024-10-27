function output_values = iiBoundaries(inputs, varargin)
% *function outputs = iiBoundaries(inputs, outputs, opts)*
%
% iiBoundaries function computes the intersection information and the angle
% (theta) between the decision boundaries of stimuli and choices based on neural
% responses.
%
% Inputs:
%   - inputs: A cell array containing the input time series data. Each cell represents a time series, where:
%             - inputs{1}: First time series (S) with dimensions        | Must be binary Stimulus
%                          nDims X nTimepoints X nTrials 
%             - inputs{2}: Second time series (R) with dimensions
%                          nDims X nTimepoints X nTrials 
%             - inputs{3}: Third time series (C) with dimensions        | Must be binary Choice
%                          nDims X nTimepoints X nTrials 
%
%   - outputs: A cell array of strings specifying the requested outputs:
%              - 'II'          : Intersection information between S R and C.
%              - 'theta'       : Angle (in degrees) between the decision boundaries of S and C.
%              - 'labelsStim'  : Predicted labels for Stimulus S by the SVM based on neural responses.
%              - 'labelsChoice': Predicted labels for Choice C by the SVM based on neural responses
%
%   - varargin: Optional arguments, passed as a structure. Fields may include:
%              - svm_options:        Specify the options that are used for the SVM 
%                                    (type 'help svm_wrapper' for more information)   
%                                    
%              - ii_options:         Specify the options that are used to compute the 
%                                    intersection information 
%                                    (type 'help II' for more information) 

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

if nargin < 1
    error('iiBoundaries:notEnoughInput', 'Not enough input arguments. See help iiBoundaries for usage info.');
elseif nargin == 1
    outputs = {'all'};
    opts = struct();
elseif nargin == 2
    if iscell(varargin{1})
        outputs = varargin{1};
        opts = struct();
    elseif isstruct(varargin{1})
        opts = varargin{1};
        outputs = {'all'};
    end
elseif nargin == 3
    outputs = varargin{1};
    opts = varargin{2};
end

if length(inputs) ~= 3
    error('iiBoundaries:InvalidInput', 'Invalid Input. See help iiBoundaries for usage info');
end

S = inputs{1};
R = inputs{2};
C = inputs{3};

if numel(unique(S)) ~= 2
    error('Stimulus must be binary to compute iiBoundaries.');
end
if numel(unique(C)) ~= 2
    error('Choice must be binary to compute iiBoundaries.');
end

nNeurons = size(R, 1);

outputs_svm = {'labels', 'mean_betaWeigths', 'intercept'};
opts.svm_options.svm_family = 'linear';
stim_SVM = svm_wrapper({R, S}, outputs_svm, opts.svm_options);
choice_SVM = svm_wrapper({R, C}, outputs_svm, opts.svm_options);

dec_SC = reduce_dim(cat(1, stim_SVM{1}, choice_SVM{1}));
II_value = II({stim_SVM{1}, dec_SC, choice_SVM{1}}, {'II(A,B,C)'}, opts.ii_options);

if nNeurons == 2
    slope_s = mean(-stim_SVM{2}(1) ./ stim_SVM{2}(2));
    slope_c = mean(-choice_SVM{2}(1) ./ choice_SVM{2}(2));
    theta = atan(abs((slope_s - slope_c) / (1 + slope_c * slope_s)));
else
    normal_vector_s = mean(stim_SVM{2});
    normal_vector_c = mean(choice_SVM{2});
    cos_theta = dot(normal_vector_s, normal_vector_c) / (norm(normal_vector_s) * norm(normal_vector_c));
    theta = acos(cos_theta);
end

theta = rad2deg(theta);

possibleOutputs = {'II', 'theta', 'labelsStim', 'labelsChoice'};
if ismember('all', outputs)
    outputs = possibleOutputs;
end

[isMember, indices] = ismember(outputs, possibleOutputs);
if any(~isMember)
    nonMembers = outputs(~isMember);
    error('iiBoundaries:invalidOutput', 'Invalid Outputs: %s', strjoin(nonMembers, ', '));
end

output_values = cell(1, length(outputs));
for i = 1:length(indices)
    idx = indices(i);
    switch possibleOutputs{idx}
        case 'II'
            output_values{i} = II_value{1};
        case 'theta'
            output_values{i} = theta;
        case 'labelsStim'
            output_values{i} = stim_SVM{1}; 
        case 'labelsChoice'
            output_values{i} = choice_SVM{1};
    end
end

end
