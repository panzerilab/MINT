function dxi = dxi_KB (K, B)
  %% Usage:  dxi = dxi_KB (K, B)
  %%
  %% Calucalates the derivative of the a priori expected value of 
  %% the entropy (xi) in the prior of the form:
  %%    P[Q]~ dQ/Q Q^{B/K}
  %% with resect to B.
  %%
  %% In:
  %%   K -- cardinality of the distribution support set;
  %%   B -- B/K extra samples is added to each possible outcome, must
  %%        be the same dimensionality object as K.
  %% Out:
  %%   dxi -- derivative of the expected entropy; 
  %%        same dimensionality as K.
  %%
  %%
  %% Depends on:
  %%   polygamma.m
  %%
  %% (c) Ilya Nemenman, 2002--2006
  %% Distributed under GPL, version 2
  %% Copyright (2006). The Regents of the University of California.

  %% This material was produced under U.S. Government contract
  %% W-7405-ENG-36 for Los Alamos National Laboratory, which is
  %% operated by the University of California for the U.S.
  %% Department of Energy. The U.S. Government has rights to use,
  %% reproduce, and distribute this software.  NEITHER THE
  %% GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS
  %% OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
  %% SOFTWARE.  If software is modified to produce derivative works,
  %% such modified software should be clearly marked, so as not to
  %% confuse it with the version available from LANL.
  %%
  %% Additionally, this program is free software; you can redistribute
  %% it and/or modify it under the terms of the GNU General Public
  %% License as published by the Free Software Foundation; either
  %% version 2 of the License, or (at your option) any later version.
  %% Accordingly, this program is distributed in the hope that it will
  %% be useful, but WITHOUT ANY WARRANTY; without even the implied
  %% warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  %% See the GNU General Public License for more details.


  if ((all(size(K)==size(B))) | (length(K)==1) | (length(B)==1)) 
    dxi=polygamma(1, B+1) - 1./K.*polygamma(1, 1+B./K);
  else
    error('Dimensions of K and B mismatch.');
  end;


