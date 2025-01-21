function f = meanS2 (B,kx,nx,K)
  %% Usage:  f = meanS2 (B,kx,nx,K)
  %%
  %% Calucalates the expected (mean) value of the entropy squared
  %% in the a priori measure of the form  [d log(Q)] prod Q^beta
  %% with the data given by nx, kx. The distributions  involved are 
  %% assumed to have K possible outcomes. Note that even though 
  %% the integration variable is xi, the argument of this function 
  %% is B.
  %%
  %% In:
  %%     B   - any matrix B=K*beta, indexing a priori expectation of S
  %%           (integration variable);  
  %%     kx, nx - row vectors; exactly kx(i) bins had nx(i) counts;
  %%           IMPORTANT: bins with 0 counts are not indexed
  %%           (that is, there is no nx(i) == 0), but only for
  %%           consistency with other routines. In principle, this
  %%           routine does not need this check.
  %%     K   - number of bins; it can be Inf if the number of bins is
  %%           unknown and assumed infinite; (scalar)
  %% Out:
  %%     f   - value of the expected square of the entropy; 
  %%           same size as B.
  %%
  %% Depends on:
  %%     polygamma.m, psi.m
  %%
  %% Addition Reference: 
  %%     Wolpert and Wolf, 1995
  %%
  %% (c) Ilya Nemenman, 2002--2011 
  %% (c) Fernando Montani, 2005 (Octave to MatLab conversion)
  %% Distributed under GPL, version 2
  %% Copyright (2006). The Regents of the University of California.

  % This material was produced in part under U.S. Government contract
  % W-7405-ENG-36 for Los Alamos National Laboratory, which is
  % operated by the University of California for the U.S.
  % Department of Energy. The U.S. Government has rights to use,
  % reproduce, and distribute this software.  NEITHER THE
  % GOVERNMENT NOR THE UNIVERSITY MAKES ANY WARRANTY, EXPRESS
  % OR IMPLIED, OR ASSUMES ANY LIABILITY FOR THE USE OF THIS
  % SOFTWARE.  If software is modified to produce derivative works,
  % such modified software should be clearly marked, so as not to
  % confuse it with the version available from LANL.
  %
  % Additionally, this program is free software; you can redistribute
  % it and/or modify it under the terms of the GNU General Public
  % License as published by the Free Software Foundation; either
  % version 2 of the License, or (at your option) any later version.
  % Accordingly, this program is distributed in the hope that it will
  % be useful, but WITHOUT ANY WARRANTY; without even the implied
  % warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
  % See the GNU General Public License for more details.

  % we require that none nx(i)=0 so that there's no problems with
  % overflow in the case of too many bins with zero count 

  if any(nx<=0) error('Bins with zero count encountered.'); end; 

  % other sanity checks
  N  = sum (kx.*nx);		% total number of samples
  if (N<=1) warning(['Too few data samples: N=' int2str(N)]); end;
  if (length(K)~=1) error('Wrong size of variable K.'); end; 

  K1 = sum (kx);		    % number of bins with nozero occupancy

  [row, col] = size(B);
  B = B(:);
  f = zeros(size(B));

  ki = kx';

  for ii=1:row*col
    Bi=B(ii);
    b = Bi/K;
    nxb = nx+b;
    nib = nxb';
    pnxb1 = psi(nxb+1);
    pnib1 = pnxb1';
    
    p0NB2 = psi(N+Bi+2);
    p1NB2 = polygamma(1, N+Bi+2);

    %------------------------------------------------------
    % i, j term, summing over all i and j (including i==j terms,
    % thus overcounting, and then correcting for it) 
    
    % ni*nj ~= 0 contribution
    f(ii) = sum(sum((nib.*(pnib1-p0NB2).*ki) * ...
                (nxb.*(pnxb1-p0NB2).*kx) - ...
                (nib.*ki)*(nxb.*kx)*p1NB2));

    % ni*b contribution    
    f(ii) = f(ii) + 2*Bi*(1-K1/K)*nxb.*((pnxb1-p0NB2)*(psi(b+1)-p0NB2) - p1NB2)*ki;
    
    % b*b contribution 
    f(ii) = f(ii) + (1-K1/K)*(1-(K1+1)/K)*Bi^2*((psi(b+1)-p0NB2)^2-p1NB2);

    %correcting for the overcounting
    f(ii) = f(ii) -(nxb.*(pnxb1-p0NB2)).^2*ki + nxb.*nxb*ki*p1NB2;
    
    %-----------------------------------------------------
    % i term
    
    % ni contribution
    f(ii) = f(ii) + (nxb.*(nxb+1).* ((psi(nxb+2) - p0NB2).^2 + ...
			     polygamma(1,nxb+2) - p1NB2))*ki;
   
    f(ii) = f(ii) + Bi*(1-K1/K)*(1+b) * ((psi(2+b)-p0NB2)^2 + polygamma(1,b+2) ...
					 - p1NB2);

    %-----------------------------------------------------
    % normalizing
    f(ii) = f(ii)/((N+Bi)*(N+Bi+1));

  end; 
  
  f(B==Inf) = log(K)^2;

  f=reshape(f, row, col);

 
