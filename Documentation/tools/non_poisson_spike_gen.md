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
