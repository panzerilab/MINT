%%% *[marginalized_pdf] = marginalize_pdf(pdf,survivingdims)*
%%%
%%% ### Description
%%% This function reduces the number of variables in a multivariate pdf to the ones specified in survivingdims. E.g. for a 3-variate input pdf with survivingdims = [1,2] the resulting reduced_pdf will be the joint probability (1,2).
%%%
%%% ### Inputs:
%%% - *pdf*: multi-variate input pdf, specified as a multi-dimensional array such as each component of `size(pdf)` is equal to the length of the corresponding random variable.
%%% - *survivingdims*: vector (or single value) of variables that will be used in the resulting reduced_pdf. If a single value the function returns the corresponding marginal pdf of such variable.
%%% - *iscont*: vector of boolean values of same length of `length(size(pdf))` (number of random variables in original pdf). True if i-th variable is continuous, false if it is discrete.
%%% - *varargin*: vectors of $$x_i$$ used to sample the input pdf. Used for the integration of the pdf over the reduced variables. There should be a $$x_i$$ vector for each of the dimensions of the input joint pdf.
%%%
%%% ### Outputs:
%%% - *marginalized_pdf*: output joint (or marginal if survivingdims is not a vector) probability between variables in survivingdims.
