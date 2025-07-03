function [corrected_v, plugin_v, shuff_all, addOut] = resampling_correction(inputs, outputs, corefunc, varargin)
% resampling_correction - Compute bias-corrected and plugin information-theoretic values
%
% Inputs:
%   - inputs: A cell array containing the input variables, where each cell corresponds to a different input source.
%             The data can be structured as either:
%             - nDims X nTrials
%             - nDims X nTimepoints X nTrials
%
%   - outputs: A cell array specifying the information measures to compute.
%
%   - corefunc: A function handle representing the core information-theoretic computation (e.g., @PID, @MI, @FIT, @II).
%
%   - varargin: Optional arguments as a structure (`opts`), including fields:
%              - shuff: The number of shuffling iterations for bias correction (default: 20).
%              - pid_constrained: A flag indicating whether to use a constrained PID (true/false) (default: true).
%              - chosen_atom: Specifies which PID atom ('Red', 'Syn', 'Unq1', 'Unq2') or other measure to compute (default: 'Syn').
%              - inputs_nD: Higher-dimensional input data for specific cases (if needed).
%
% Outputs:
%   - corrected_v: A cell array containing the bias-corrected values for each specified output, computed by subtracting
%                  the shuffled contributions from the plugin values.
%   - plugin_v: A cell array containing the plugin (uncorrected) values, as computed from the original input-output
%              relationships.
%   - shuff_all: A cell array containing the results of the shuffled computations across all trials and outputs,
%                used for bias subtraction.
%   - addOut: Additional outputs specific to the core function used (e.g., decomposed PID atoms, FIT interaction measures).
%
% Note:
% Users can customize the core calculation function to extend the analysis for specific types of data
% or information-theoretic measures. The function is designed to handle various dimensions of input data,
% ensuring flexibility for different experimental setups.

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

opts = varargin{1};
if nargin > 5
    inputs_nD = varargin{2};
end
addOut = 0;
plugin_opts = opts;
plugin_opts.bias = 'plugin';
plugin_opts.bin_method = {'none'};
nVars = length(inputs);
nTimepoints_all = 1;
Dims_tmp = size(inputs{1});
for var = 1:nVars
    Dims_var = size(inputs{var});
    nTrials = Dims_var(end);
    if length(Dims_var) == 3
        nTimepoints_all = [nTimepoints_all, Dims_var(2)];
    end
end
nTimepoints = max(nTimepoints_all);
corrected_v = cell(1, length(outputs));


shuff_all = repmat({zeros(plugin_opts.shuff, nTimepoints)},1, length(outputs));
plugin_v = feval(corefunc, inputs, outputs, plugin_opts);
if length(inputs)==2
    meas_pdist = prob_estimator(inputs,{'P(A,B)'},plugin_opts);
else
    meas_pdist = prob_estimator(inputs,{'P(A,B,C)'},plugin_opts);
end
noutputs = length(outputs);

if plugin_opts.parallel
    shuffall_par = zeros(noutputs, plugin_opts.shuff, nTimepoints);
    parfor sIdx = 1:plugin_opts.shuff
        inputs_sh = inputs;
        if length(size(inputs_sh{1}))==3        % Generate uniform random numbers
            for it=1:nTimepoints
                cdf = cumsum(meas_pdist{it,1}(:)); % Convert probability matrix to a cumulative distribution function (CDF)
                rand_vals = rand(nTrials, 1);
                % Find corresponding indices
                [~, sampled_indices] = histc(rand_vals, [0; cdf]);
                if length(inputs_sh)==2
                    % Convert linear indices to subscripts
                    [xx, yy] = ind2sub(size(meas_pdist{1}), sampled_indices);
                    inputs_sh{1}(:,it,:) = xx;
                    inputs_sh{2}(:,it,:) = yy;
                elseif length(inputs_sh)==3
                    [xx1, xx2, yy] = ind2sub(size(meas_pdist{1}), sampled_indices);
                    inputs_sh{1}(:,it,:) = xx1;
                    inputs_sh{2}(:,it,:) = xx2;
                    inputs_sh{3}(:,it,:) = yy;
                end
            end
        else
            cdf = cumsum(meas_pdist{1}(:)); % Convert probability matrix to a cumulative distribution function (CDF)
            rand_vals = rand(nTrials, 1);
            % Find corresponding indices
            [~, sampled_indices] = histc(rand_vals, [0; cdf]);
            if length(inputs_sh)==2
                % Convert linear indices to subscripts
                [xx, yy] = ind2sub(size(meas_pdist{1}), sampled_indices);
                inputs_sh{1}(:,:) = xx;
                inputs_sh{2}(:,:) = yy;
            elseif length(inputs_sh)==3
                [xx1, xx2, yy] = ind2sub(size(meas_pdist{1}), sampled_indices);
                inputs_sh{1}(:,:) = xx1;
                inputs_sh{2}(:,:) = xx2;
                inputs_sh{3}(:,:) = yy;
            end
        end    
        shuff_v  = feval(corefunc, inputs_sh, outputs, plugin_opts);
        for outIdx = 1:noutputs
            shuffall_par(outIdx,sIdx,:) = shuff_v{outIdx};
        end
    end
    for spar=1:noutputs
        shuff_all{spar} = reshape(shuffall_par(spar,:,:),plugin_opts.shuff, nTimepoints);
    end
else
    for sIdx = 1:plugin_opts.shuff
        inputs_sh = inputs;
        if nTimepoints>1        % Generate uniform random numbers
            for it=1:nTimepoints
                cdf = cumsum(meas_pdist{it,1}(:)); % Convert probability matrix to a cumulative distribution function (CDF)
                rand_vals = rand(nTrials, 1);
                % Find corresponding indices
                [~, sampled_indices] = histc(rand_vals, [0; cdf]);
                if length(inputs_sh)==2
                    % Convert linear indices to subscripts
                    [xx, yy] = ind2sub(size(meas_pdist{it,1}), sampled_indices);
                    inputs_sh{1}(:,it,:) = xx;
                    inputs_sh{2}(:,it,:) = yy;
                elseif length(inputs_sh)==3
                    [xx1, xx2, yy] = ind2sub(size(meas_pdist{it,1}), sampled_indices);
                    inputs_sh{1}(:,it,:) = xx1;
                    inputs_sh{2}(:,it,:) = xx2;
                    inputs_sh{3}(:,it,:) = yy;
                end
            end
        else
            cdf = cumsum(meas_pdist{1}(:)); % Convert probability matrix to a cumulative distribution function (CDF)
            rand_vals = rand(nTrials, 1);
            % Find corresponding indices
            [~, sampled_indices] = histc(rand_vals, [0; cdf]);
            if length(inputs_sh)==2
                % Convert linear indices to subscripts
                [xx, yy] = ind2sub(size(meas_pdist{1}), sampled_indices);
                inputs_sh{1} = xx';
                inputs_sh{2} = yy';
            elseif length(inputs_sh)==3
                [xx1, xx2, yy] = ind2sub(size(meas_pdist{1}), sampled_indices);
                inputs_sh{1} = xx1';
                inputs_sh{2} = xx2';
                inputs_sh{3} = yy';
            end
        end
        
        shuff_v  = feval(corefunc, inputs_sh, outputs, plugin_opts);
        for outIdx = 1:noutputs
            shuffall_par(outIdx,sIdx,:) = shuff_v{outIdx};
        end
    end
    for spar=1:noutputs
        shuff_all{spar} = reshape(shuffall_par(spar,:,:),plugin_opts.shuff, nTimepoints);
    end
end
for outIdx = 1:length(outputs)
    corrected_v{outIdx} = 2*plugin_v{outIdx} - mean(shuff_all{outIdx},1);
end
end
