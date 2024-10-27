%%% *function [p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] = build_p(S, R, C)*
%%%
%%% ### Description
%%% The function estimates the probability distribution p(s,r,c) from the 3D histogram of the input (S, R, C) occurrences.
%%%
%%% ### Inputs:
%%% - *S*: must be a one dimensional array of *n_trials X 1* elements representing the discrete value of the stimulus presented in each trial.
%%% - *R*: must be a one dimensional array of *n_trials X 1* elements representing the response at each trial, here we suppose the response has been binned already.
%%% - *C*: must be a one dimensional array of *n_trials X 1* elements representing the discrete value of the choice made by the subject in each trial.
%%%
%%% ### Outputs:
%%% - *p_src*: joint probability p(s,r,c).
%%% - *p_crs*: joint probability p(c,r,s).
%%% - *n_S*: number of stimuli.
%%% - *n_R*: number of responses.
%%% - *n_C*: number of choices.
%%% - *n_singleton_dims*: number of singleton dimensions.
