%%% *function [pid_v, table_prob] = pidBROJA(pdf_dirty)*
%%%
%%% ### Description
%%% Compute the Partial Information Decomposition (PID) according to the BROJA-2 method. The function evaluates the unique, shared,
%%% and synergistic information contributions between two input variables (X, Y) with respect to a target variable (Z).
%%%
%%% ### Inputs:
%%% - pdf_dirty: A probability distribution function array (nX x nY x nZ), representing the joint probability of the variables X, Y, and Z
%%%
%%% ### Outputs:
%%% - pid_v: A vector containing the information components in the following order: Shared information (SI), unique information of Y (UIY), unique information of Z (UIZ), and complementary information (CI).
%%% - table_prob: The non-zero probabilities used in the computation
%%%
%%% ### Further notes:
%%% The function uses the BROJA-2 method to optimize the PID measures.
