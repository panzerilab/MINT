function [spikes] = poisson_spike_gen(time, rate, noise_prob)
%%% *function [spikes] = poisson_spike_gen(time, rate, noise_prob)*
%%%
%%% ### Description
%%% poisson_spike_gen generates a train of Poisson spikes with defined rate.
%%%
%%% ### Inputs:
%%% - *time*: time array for the spike generation of size *n_timesteps*.
%%% - *rate*: rate for the Poisson process, if a double it is assumed a constant rate, if a vector a time varying rate can be specified. In the latter case *length(rate) = n_timesteps*.
%%% - *noise_prob*: noise probability of Poisson process (probability of a spike to be randomly generated/suppressed.
%%%
%%% ### Outputs:
%%% - *spikes*: spike train corresponding to given inputs.

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

if length(rate)==1
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

