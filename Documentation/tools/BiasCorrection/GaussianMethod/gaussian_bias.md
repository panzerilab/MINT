%%% *function bias = gaussian_bias(Nt,L)*
%%%
%%% ### Description
%%% Computes bias of gaussian entropy estimates.
%%%
%%% ### Inputs:
%%% - *Nt*: number of trials used for the estimation of entropy. Can be a integer scalar or a column array of integers: in this case the bias is computed for each of the values in Nt.
%%% - *L*: dimensionality of the response.
%%%
%%% ### Outputs:
%%% - *bias*: bias of gaussian entropy estimates.
%%%
