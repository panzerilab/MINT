%%% *function pars = build_parameters_structure(X, opts, responseMatrixName, outputsList)*
%%%
%%% ### Description
%%% Build parameters structure for ENTROPY.m
%%%
%%% ### Inputs:
%%% - *X*: *nDimensionsX X nTrials X nY* 3D array including response data.
%%% - *opts*: options to compute the entropy.
%%% - *responseMatrixName*: name associated to the X matrix.
%%% - *outputsList*: cell array of char arrays of strings specifying which quantities to compute.
%%%
%%% ### Further notes:
%%%
%%% Recap of the fields in the parameter structure
%%% ----------------------------------------------
%%%
%%% X-related parameters:
%%%   - pars.Nc
%%%   - pars.Ns
%%%   - pars.Nt
%%%   - pars.maxNt
%%%   - pars.totNt
%%%
%%% TESTMODE-related parameters:
%%%   - pars.testmode
%%%
%%% METHOD-related parameters:
%%%   - pars.methodFunc
%%%   - pars.methodNum
%%%
%%% EXTRAPOLATION-related parameters
%%%   - pars.xtrp
%%%
%%% BIAS-related parameters:
%%%   - pars.biasCorrNum
%%%   - pars.biasCorrFuncName
%%%
%%% BOOTSTRAP-related parameters:
%%%   - pars.btsp
%%%
%%% OUTPUT-related parameters:
%%%   - pars.whereHX
%%%   - pars.whereHXY
%%%   - pars.whereHlX
%%%   - pars.whereHlXY
%%%   - pars.whereHiX
%%%   - pars.whereChiX
%%%   - pars.whereHshX
%%%   - pars.whereHshXY
%%%   - pars.whereHiXY
%%%   - pars.doHX
%%%   - pars.doHXY
%%%   - pars.doHlX
%%%   - pars.doHlXY
%%%   - pars.doHiX
%%%   - pars.doHiXY
%%%   - pars.doChiX
%%%   - pars.doHshX
%%%   - pars.doHshXY
%%%
%%% CHECKS-related parameters:
%%%   - pars.addChecks
%%%   - pars.numberOfSpecifiedOptions
%%%   - pars.Noutput
%%%
