function f = meanS (B,kx,nx,K)
  %% Usage:  f = meanS (B,kx,nx,K)
  %%
  %% Calucalates the expected (mean) value of the entropy 
  %% in the a priori measure of the form  [d log(Q)] ...prod Q^beta
  %% with the data given by nx, kx. The distributions  involved are 
  %% assumed to have K possible outcomes. Note that even though 
  %% the integration variable is ...xi, the argument of this function 
  %% is B.
  %%
  %% In:
  %%     B   - any matrix of B=K*...beta, indexing a priori expectation of S
  %%           (integration variable);  
  %%     kx, nx - row vectors; exactly kx(i) bins had nx(i) counts;
  %%           IMPORTANT: bins with 0 counts are not indexed
  %%           (that is, there is no nx(i) == 0), but only for
  %%           consistency with other routines. In principle, this
  %%           routine does not need this check.
  %%     K   - number of bins; it can be Inf if the number of bins is
  %%           unknown and assumed infinite; (scalar)
  %% Out:
  %%     f   - value of the expected entropy; same size as B.
  %%
  %% Depends on:
  %%     B_xiK.m, psi.m
  %%
  %% Addition Reference: 
  %%     Wolpert and Wolf, 1995
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

  %% we require that none nx(i)=0 so that there's no problems with
  %% overflow in the case of too many bins with zero count 

  if any(nx<=0) error('Bins with zero count encountered.'); end;

  %% other sanity checks
  N  = sum (kx.*nx);		% total number of samples
  if (N<=1) warning(['Too few data samples: N=' int2str(N)]); end;
  if (length(K)~=1) error('Wrong size of variable K.'); end; 

  K1 = sum (kx);		% number of bins with nozero occupancy

  [row, col] = size(B);
  B = B(:);


  ovrNB = 1./(N+B);
  osB = ones(size(B));
  osn = ones(size(nx));

  f=zeros(size(B));
  %% sum is written in a form to avoid the lossof precision as K->Inf
  f = psi(N+B+1) - (ovrNB*osn).*(osB*nx + ...
			  B/K*osn).*psi(osB*nx+B/K*osn+1)*kx(:) - ...
      B.*ovrNB*(1-K1/K).*psi(1+B/K);

  f(find(B==Inf)) = log(K);


  f=reshape(f,row,col);


 
  
