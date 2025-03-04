function [spikes] = non_poisson_spike_gen(time, rate, noise_prob)
% non_poisson_spike_gen - Generate non-Poisson spike trains using a gamma-distributed inter-spike interval (ISI)
%
% This function generates a spike train where the inter-spike intervals (ISI) follow a gamma distribution. The
% function can generate spike trains with different levels of regularity by controlling the shape parameter of 
% the gamma distribution. It allows generating super-Poisson, Poisson-like, or sub-Poisson spike patterns.
%
% Inputs:
%   - time: A vector specifying the time points for spike generation. Must be of size *n_timesteps*.
%           It defines the temporal resolution and the duration of the simulation.
%   - rate: A scalar value specifying the firing rate of the spike train (in spikes per unit time).
%           This defines the expected rate of spike generation.
%   - noise_prob: A scalar value representing the shape parameter of the gamma distribution:
%                 - If *noise_prob* < 1, the generated spikes will exhibit super-Poisson variability (more irregular).
%                 - If *noise_prob* = 1, the spikes will follow a Poisson process (random, no noise correlation).
%                 - If *noise_prob* > 1, the spikes will exhibit sub-Poisson variability (more regular).
%
% Outputs:
%   - spikes: A binary vector (1 x *n_timesteps*) representing the generated spike train. Each element corresponds
%             to a time step in the input *time* array, where a 1 indicates a spike at that time step, and a 0 indicates
%             no spike.
%
% Method:
% The function models the inter-spike intervals (ISI) using a gamma distribution with shape parameter *noise_prob* 
% and a mean ISI equal to the inverse of the input *rate*. The spike times are calculated by summing successive ISIs,
% and spikes are placed at the closest time points in the *time* array. The output is a spike train vector that can
% be used in neural simulations or analyses.
%
% Example:
% Suppose you want to generate a spike train over 10 seconds with a firing rate of 5 Hz:
%
%     time = 0:0.001:10;  % Time vector with 1 ms resolution
%     rate = 5;           % Firing rate in Hz
%     noise_prob = 2;     % shape parameter
%
%     spikes = non_poisson_spike_gen(time, rate, noise_prob);
%
% This will return a spike train of size (1 x 10001) corresponding to the time array. The spike intervals will be
% gamma-distributed with a shape parameter of 2, resulting in more regular, less noisy spikes compared to a Poisson process.

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
    error('nonpoissonspikegen:notEnoughInput', msg);
end
if any(isnan(time))
    msg = "time contains NaNs. Aborting.";
    error('nonpoissonspikegen:NaNInput', msg);
end
if isnan(rate)
    msg = "rate contains NaNs. Aborting.";
    error('nonpoissonspikegen:NaNInput', msg);
end
if isnan(noise_prob)
    msg = "noise_prob contains NaNs. Aborting.";
    error('nonpoissonspikegen:NaNInput', msg);
end

dt = diff(time); dt = dt(1);
assert(isscalar(rate));
spikes = zeros(1,length(time));

last_spike_time = 0;
mean_isi = 1/rate;
theta = mean_isi/noise_prob;
spike_times = [];
while last_spike_time < time(end)
    isi = gamrnd(noise_prob,theta);
    last_spike_time = last_spike_time + isi;
    if last_spike_time <= time(end)
        spike_times = [spike_times last_spike_time];
    end
end
% create spike train
for j=1:length(spike_times)
    spikes(round(spike_times(j)/dt)+1) = 1;
end

end

