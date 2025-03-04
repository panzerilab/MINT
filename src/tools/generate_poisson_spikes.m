function spikeTrains = generate_poisson_spikes(nNeurons,tPoints,nTrials_per_stim,lambda, varargin)
% function spikeTrains = generate_poisson_spikes(nNeurons, tPoints, nTrials_per_stim, lambda, varargin)
%
% The generate_poisson_spikes function generates spike train data for a population of neurons 
% using Poisson processes, incorporating stimulus-specific firing rates and optional background noise.
% This is useful for simulating neural activity in response to different stimuli for computational 
% neuroscience analyses.
%
% Inputs:
%   - nNeurons: Integer specifying the number of neurons to simulate.
%
%   - tPoints: Integer representing the number of time points (or time bins) in each trial. 
%              This determines the temporal resolution of the spike trains.
%
%   - nTrials_per_stim: Integer defining the number of trials for each stimulus. 
%                       The function generates multiple trials for statistical reliability.
%
%   - lambda: A vector of length nStimuli, specifying the mean firing rate (in Hz) for each stimulus 
%             condition. Each element corresponds to the firing rate for one stimulus.
%
%   - varargin: Optional arguments passed as additional parameters. 
%               - lambda_noise (optional): A vector of length nStimuli, specifying the background 
%                 noise firing rate (in Hz) for each stimulus. If not provided, the function assumes 
%                 no noise (default: zeros(1, nStimuli)).
%
% Outputs:
%   - spikeTrains: A 3D matrix of size [nNeurons x tPoints x (nTrials_per_stim * nStimuli)], where:
%                  - The first dimension corresponds to neurons.
%                  - The second dimension corresponds to time points.
%                  - The third dimension contains spike data for all trials and stimuli, concatenated.
%
%                  The entries in the matrix are binary (0 or 1), where 1 indicates a spike at a 
%                  particular time point for a given neuron and trial.
%
% Example:
% Generate spike trains for 10 neurons, over 100 time points, with 5 trials per stimulus, 
% for 2 stimuli with firing rates [20, 50] Hz, and background noise rates [5, 10] Hz:
% nNeurons = 10;
% tPoints = 100;
% nTrials_per_stim = 5;
% lambda = [20, 50]; % Stimulus-driven firing rates
% lambda_noise = [5, 10]; % Background noise rates
% spikeTrains = generate_poisson_spikes(nNeurons, tPoints, nTrials_per_stim, lambda, lambda_noise);
%
% Limitations:
% - Poisson processes assume independence of spike occurrences, which may not always hold for real neural data.

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

nStimuli = length(lambda);
spikeTrains = zeros(nNeurons, tPoints, nTrials_per_stim*nStimuli);
if ~isempty(varargin) && length(varargin) >= 1
        lambda_noise = varargin{1}; 
    else
        lambda_noise = zeros(1, nStimuli); 
    end
for stim = 1:nStimuli
    R = zeros(nNeurons,tPoints,nTrials_per_stim);
    lambda_stim = lambda(stim);
    lambda_noise_stim = lambda_noise(stim);
    for trialIdx = 1:nTrials_per_stim
        spikes_noise = poisson_spike_gen(1:tPoints, lambda_noise_stim/tPoints, 0);
        for cellIdx = 1:nNeurons
            spikes = poisson_spike_gen(1:tPoints, lambda_stim/tPoints, 0);
            R(cellIdx,:,trialIdx) = spikes + spikes_noise;           
        end
    end
    startIdx = (stim - 1) * nTrials_per_stim + 1;
    endIdx = stim * nTrials_per_stim;
    spikeTrains(:, :, startIdx:endIdx) = R;
end


