%%% *function ind = findtrial(Nt, maxNt, Ns)*
%%%
%%% ### Description
%%% Index responses corresponding to trials.
%%%
%%% ### Inputs:
%%% - *Nt*: number of trials per stimulus array.
%%% - *maxNt*: max number of trials (must be equal to max(NT)).
%%% - *Ns*: number of stimuli (must be equal to length(NT)).
%%%
%%% ### Outputs:
%%% - *ind*: position of the trials in the response matrix R (see ENTROPY and INFORMATION) ignoring the first dimension (dimensionality of the response).
%%%
%%% ### Further notes:
%%%
%%% EXAMPLE
%%%   Let's consider the response matrix in the figure below (letters are
%%%   used to represent undefined values):
%%%
%%%                             ----------------- 
%%%                             | q | s | - | - |
%%%                             |---------------|
%%%                             | r | t | - | - |
%%%                             -----------------
%%%                        -----------------        /
%%%                        | i | m | o | - |       /
%%%                        |---------------|      /
%%%                        | l | n | p | - |     / 
%%%                        -----------------    /
%%%              |    -----------------        / S=3
%%%              |    | a | c | e | g |       /
%%%          L=2 |    |---------------|      /
%%%              |    | b | d | f | h |     /
%%%              |    -----------------    /
%%%      
%%%                   -----------------
%%%                          T=4
%%%
%%%   in this case L=2, T=4 and S=3. Four trials are availabe for the first
%%%   stimulus, two for the second and three for the third stimulus. We have
%%%   Nt = [4; 3; 2] and maxNt = 4 and thus
%%%
%%%         ind = tridxs([4; 3; 2], 4, 3)
%%% 
%%%         ind =
%%% 
%%%              1
%%%              2
%%%              3
%%%              4
%%%              5
%%%              6
%%%              7
%%%              9
%%%             10
%%%
%%%   The indexes computed this way can be used to extract all the responses
%%%   corresponding to trials by means of the command
%%%         
%%%          R(:,ind) = 
%%%
%%%              a      c      e      g      i      m      o      q      s
%%%              b      d      f      h      l      n      p      r      t
%%%
%%% REMARKS
%%%   The inputs are obviously redundant since all the information required for
%%%   the computation is contained in NT. However, by passing the additional
%%%   arguments we can speed up things a bit.
%%%
%%%
