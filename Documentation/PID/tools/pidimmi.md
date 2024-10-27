%%% *function pid_v = pidimmi(pdf_dirty)*
%%%
%%% ### Description
%%% Compute Partial Information Decomposition (PID) values from a given probability distribution (pdf_dirty).
%%%
%%% Partial Information Decomposition separates the information that multiple sources provide about a target into unique, redundant, and synergistic components. 
%%% This function takes a 3D probability distribution (joint pdf of three variables) and calculates the PID values.
%%%
%%% ### Inputs:
%%% - pdf_dirty: A 3-dimensional array representing the joint probability density function of three variables, typically denoted as X, Y, and Z. 
%%%
%%% ### Outputs:
%%% - pid_v: A vector containing four values representing the components of PID:
%%% - Redundant information (common information both sources provide about the target)
%%% - Unique information from the first variable (information only the first source provides about the target)
%%% - Unique information from the second variable (information only the second source provides about the target)
%%% - Synergistic information (information that is only available when both sources are considered together)
