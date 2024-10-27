%%% *function windowVec = windowing(time,nWindows)*
%%% 
%%% ### Description
%%% windowing chooses the best way of splitting the data given the amount of windows provided by the user. The amount of windows can be slightly modified from the input.
%%%
%%% ### Inputs:
%%% - *time*: vector including time steps.
%%% - *nWindows*: integer including how many windows we want to split the data into.
%%%
%%% ### Outputs:
%%% - *windowVec*: vector inluding the length of each window.
