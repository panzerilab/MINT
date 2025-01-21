function f = mlog_evidence(B, kx, nx, K)
  %% Usage:  f = mlog_evidence (B, kx, nx, K)
  %%
  %% Computes the 'action' (the negative logarithms of the 
  %% evidence) for the integral over ...xi (see NSB method for 
  %% calculating entropies of discrete pdf's). Does not include
  %% the contribution from the prior over ...xi. Note that even though 
  %% the integration variable is ...xi, the argument of this function 
  %% is B.
  %%
  %% In:
  %%     B   - any matrix, the a priori expectation of S
  %%           (integration variable);  
  %%     kx, nx - row vectors; exactly kx(i) bins had nx(i) counts;
  %%           IMPORTANT: bins with 0 counts are not indexed
  %%           (that is, there is no nx(i) == 0) 
  %%     K   - number of bins; it can be Inf if the number of bins is
  %%           unknown and assumed infinite; 
  %% Out:
  %%     f   - value of the negative log-evidence(action); same 
  %%           size as B.
  %%
  %% Depends on:
  %%     B_xiK.m, lpoch.m
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
  if any(nx<=0) error('Bins with zero count encountered'); end;

  %% other sanity checks
  N  = sum (kx.*nx);		% total number of samples
  if (N<=1) error(['Too few data samples: N=' int2str(N)]); end; 
  if (length(K)~=1) error('Wrong size of variable K.'); end; 

  [row, col] = size(B);
  B = B(:);
  f = zeros(size(B));
  B0 =find(B==0);		% B's that are 0

  %% if the are zeros in B, call mlog_evidence recursively with
  %% only those B's that are not zero, and state f(B==0)=Inf.
  if (length(B0)>0)
    Bn0=find(B~=0);		% B's that are not 0
    f(Bn0) = mlog_evidence(B(Bn0), kx, nx, K);
    f(B0) = Inf;
  else
    %% else just evaluate f(B) 
    ng1= find(nx>1);		% occupancy which is greater than 1
    K1 = sum (kx);		% number of bins with nozero occupancy
    K2 = sum (kx(ng1));		% number of bins with coincidences

    if (all(round(nx)==nx))	% all integer occupancies
      ints = 1;
    else			% fractional occupancies
      ints = 0;
    end

    if(ints==1)			% integer counts
      if(length(ng1)>0)		% if there are coincidences
	f1 =  K2*gammaln(1+B/K) ;
    f2 = -gammaln(B*ones(1,length(ng1))/K + ones(length(B),1)*nx(ng1));
    f3=  f2*(kx(ng1))';

    f = f1 + f3(:);
    

   
      end;
    else			
     
      f = -(gammaln(B*ones(1,length(nx))/K + ...
		      ones(length(B),1)*nx)*kx(:)) ...
	    + K1*gammaln(1+B/K);
    end;
 
     
    %% Need to calcuate f += -K1*log(B) +gammaln(B+N) - gammaln(B)
    %% But to avoid (big number) - (big number) = lost precision
    %% problem, need to treat different regimes of N and B
    %% differently
    
    %% First regime, aymptotically large B and B/N: B>100 and N<0.01*B;
    %% here we can expand gammaln-gammaln = psi*N + psi_1/2*N^2+...
    large = (B> max([100, 100*N]));


    %% polygamma(n,B) ~B^(-n); thuskx(ng1) we expand in (N/B);
    %% which of expansion parameters is the worst? how many 
    %% series terms will we need? as seen below, the leading
    %% term in f for N==K1 (worst case) is psi_asymp(B)/N ~ N/B.
    %% we need 10^(-15) precision relative to that term. Note
    %% that the series expansion has the form 
    %% f = leading + psi_1/2!*N^2 +psi_2/3!*N^3 +... =
    %%   = leading + (N/B +(N/B)^2+...)
    if(sum(large)>0)
      nterms  = max(ceil(abs((-15 -log10(N))./log10(N./B(large))))) + 1;
      ifac=1./cumprod([2:(nterms+1)]); % factorial coefficients
      %% we add starting from the smallest term
      for i = nterms:-1:1
	f(large) = f(large) + polygamma(i,B(large)).*ifac(i).*N^(i+1);      
      end; 
      f(large) = f(large) + (N-K1)*log(B(large)) + psi_asymp(B(large))*N;
    end;

    
    other = ~large;
    
    f(other) =  f(other) - K1*log(B(other)) + gammaln(B(other)+N) -gammaln(B(other));
  end;


