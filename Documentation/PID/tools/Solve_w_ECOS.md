%%% BROJA_2PID.py -- Python module
%%%
%%% BROJA_2PID: Bertschinger-Rauh-Olbrich-Jost-Ay (BROJA) bivariate Partial Information Decomposition
%%% https://github.com/Abzinger/BROJA_2PID
%%% (c) Abdullah Makkeh, Dirk Oliver Theis
%%% Permission to use and modify with proper attribution
%%% (Apache License version 2.0)
%%%
%%% Information about the algorithm, documentation, and examples are here:
%%% @Article{makkeh-theis-vicente:pidOpt:2017,
%%%          author =       {Makkeh, Abdullah and Theis, Dirk Oliver and Vicente, Raul},
%%%          title =        {BROJA-2PID: A cone programming based Partial Information Decomposition estimator},
%%%          journal =      {jo},
%%%          year =         2017,
%%%          key =       {key},
%%%          volume =    {vol},
%%%          number =    {nr},
%%%          pages =     {1--2}
%%% }
%%% Please cite this paper when you use this software (cf. README.md)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% ECOS's exp cone: (r,p,q)   w/   q>0  &  exp(r/q) <= p/q
%%% Translation:     (0,1,2)   w/   2>0  &  0/2      <= ln(1/2)
%%%
%%% Methods:
%%%   - Solve_w_ECOS: Constructor method
%%%   - create_model: Creates the optimization model
%%%   - r_vidx: Computes the index for variable r_i
%%%   - p_vidx: Computes the index for variable p_i
%%%   - q_vidx: Computes the index for variable q_i
%%%   - solve: Solves the optimization problem
%%%   - entropy_X: Computes the entropy of variable X
%%%   - entropy_X2: Computes the entropy of variable X using precomputed table
%%%   - condenttropy: Computes conditional entropy
%%%   - addqdistr: Adds quantum distribution to the triplets
%%%   - condentropy2: Computes conditional entropy using precomputed table
%%%   - condentropy_orig: Computes original conditional entropy
%%%   - create_qxyz: Creates q_xyz table
%%%   - condYmutinf: Computes conditional mutual information for Y
%%%   - condZmutinf: Computes conditional mutual information for Z
%%%   - dual_value: Computes the dual value
%%%   - check_feaseability: Checks the feasibility of the solution
