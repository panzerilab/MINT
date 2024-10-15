function [spikes] = non_poisson_spike_gen(time, rate, noise_prob)
%%% ### Description
%%% non_poisson_spike_gen generates a train of non-Poisson spikes with
%%% defined rate using a gamma distribution for the ISI
%%%
%%% ### Inputs:
%%% - *time*: time array for the spike generation of size *n_timesteps*.
%%% - *rate*: rate for the Poisson process, must be a single value.
%%% - *noise_prob*: shape parameter of the gamma distribution (noise_prob<1 super-Poisson noise_prob=1 Poisson-like noise_prob>1 sub-Poisson).
%%%
%%% ### Outputs:
%%% - *spikes*: spike train corresponding to given inputs.

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
assert(length(rate)==1);
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

