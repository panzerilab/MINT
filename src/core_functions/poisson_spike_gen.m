function [spikes] = poisson_spike_gen(time, rate, noise_prob)
% poisson_spike_gen - Generate a Poisson spike train with optional noise
%
% This function generates a Poisson-distributed spike train based on a specified firing rate. It can introduce
% noise to the spike train by randomly flipping spike states with a given probability, controlled by the *noise_prob* parameter.
% The function can handle a constant or time-varying firing rate and produces a binary spike train, where 1 indicates a spike 
% and 0 indicates no spike at the respective time step.
%
% Inputs:
%   - time: A vector specifying the time points for spike generation. The length of this vector defines the temporal
%           resolution and the duration of the spike train. Must be of size *n_timesteps*.
%   - rate: A scalar or a vector specifying the firing rate (in spikes per unit time). If a scalar is provided, the firing
%           rate is assumed to be constant over time. If a vector is provided, it should match the size of the *time* vector
%           to allow for time-varying rates.
%   - noise_prob: A scalar value between 0 and 1. This parameter controls the amount of noise introduced into the spike train
%                 by flipping spike states (1 to 0 or 0 to 1) with the specified probability:
%                 - *noise_prob = 0*: No noise is added, and the spike train follows a perfect Poisson process.
%                 - *noise_prob > 0*: The spike train is corrupted by random noise, where some spikes may be randomly
%                 added or removed.
%
% Outputs:
%   - spikes: A binary vector of size (1 x *n_timesteps*) representing the generated spike train. Each element is either 1 
%             (indicating a spike) or 0 (indicating no spike).
%
% Method:
% The function generates a Poisson spike train by comparing random numbers drawn from a uniform distribution [0, 1]
% to the probability of firing at each time step, given by *dt* (time step size) times the firing rate. If noise is introduced,
% additional random numbers are used to flip spike states with a probability equal to *noise_prob*.
%
% Example:
% Suppose you want to generate a Poisson spike train with a firing rate of 10 Hz over 5 seconds with 1 ms time resolution:
%
%     time = 0:0.001:10;  % Time vector with 1 ms resolution
%     rate = 5;           % Firing rate in Hz
%     noise_prob = 0.1;   % Introduce 10% noise
%
%     spikes = poisson_spike_gen(time, rate, noise_prob);
%
% This will return a spike train of size (1 x 10001), corresponding to the time vector, with a Poisson-distributed spike pattern
% and some random noise affecting 10% of the spike occurrences.

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

if nargin < 3
    msg = "Not enough input arguments.";
    error('poissonspikegen:notEnoughInput', msg);
end

if any(isnan(time))
    msg = "time contains NaNs. Aborting.";
    error('poissonspikegen:NaNInput', msg);
end
if isnan(rate)
    msg = "rate contains NaNs. Aborting.";
    error('poissonspikegen:NaNInput', msg);
end
if isnan(noise_prob)
    msg = "noise_prob contains NaNs. Aborting.";
    error('poissonspikegen:NaNInput', msg);
end

dt = time(2) - time(1);

if isscalar(rate)
    rate = rate*ones(size(time));
else
    assert(length(rate) == length(time))
    if size(rate) ~= size(time)
        time = time';
    end
end

spikes = zeros(size(time));
rand_nums = unifrnd(0,1,size(time));
spikes(rand_nums <= dt*rate) = 1;

if noise_prob ~= 0
    rand_nums = unifrnd(0,1,size(time));
    spikes(rand_nums < noise_prob) = ~spikes(rand_nums < noise_prob);
end

end

