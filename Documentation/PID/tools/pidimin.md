%%% function pid_v = pidimin(pdf_dirty)
%%%
%%% ### Description
%%% Compute partial information decomposition (PID) values for a given probability distribution.
%%%
%%% ### Inputs:
%%% - pdf_dirty: n X m X k array containing the joint probability distribution of three variables. The values must sum to one and cannot contain negative or NaN values.
%%%
%%% ### Outputs:
%%% - pid_v: 1x4 vector containing the PID values for different informational components:
%%%     - Synergy: shared information that only the combination of all variables provides.
%%%     - Unique1: information uniquely contributed by the first variable.
%%%     - Unique2: information uniquely contributed by the second variable.
%%%     - Redundancy: information that is redundantly contributed by any of the variables.
%%%
%%% #### Dependency
%%% This function depends on pid_lattice which calculates the PID components based on the provided probability distribution and desired information structure (like synergy, redundancy, etc.).
