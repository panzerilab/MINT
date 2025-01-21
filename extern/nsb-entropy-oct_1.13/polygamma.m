function f = polygamma (n, x)
  %% Usage:  f = polygamma (n, x)
  %% 
  %% Calculates the n'th order polygamma function of the argument x.
  %% Input:
  %%   n - a scalar, the order of the polygamma to be calculated;
  %%       0<=n<...infty, we floor n to the nearest integer;
  %%       for n=0 (digamma), psi(x) is evaluated.
  %%   x - arguments, at which polygamma is to be evaluated; any
  %%       matrix, x's may not be complex.
  %%
  %% Output:
  %%   f - value of the function, same dimensions as x
  %%
  %% Depends on:
  %%   psi.m
  %%
  %% We use the analytical derivative of the Stirling's series for the
  %% digamma function to calculate the higher order polygamma's
  %%   psi_0(x+1) = log(x) +1/2x - ...sum_{i=1}^{...infty} B_{2i}/ (2i x^{2i})
  %% that is:
  %%   psi_n(x+1) = (-1)^(n-1) (n-1)!/x^n +1/2 n! (-1)^n/x^(n+1) -
  %%              - ...sum_{i=1}^{...infty} (-1)^n 
  %%                     B_{2i} prod[2i+1:2i+n-1]/(x^{2i+n}),
  %% where B_{2i} are the corresponding Bernoulli number.
  %% 
  %% We also use the recurrence relation
  %%   psi_n(x) = psi_n(x+1) + (-1)^(n+1) n!z^(-n-1)
  %% to move the value of the argument to the asymptotic regime.
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


  %% Bernoulli numbers, b(i) is B_{2i}, up to B_40;
  %% this provides ~full machine precision
  b=[ -0.5, 0.166666666666667, 0, -0.0333333333333333,...
     0, 0.0238095238095238, 0, -0.0333333333333333,...
     0, 0.0757575757575758, 0, -0.253113553113553,...
     0, 1.16666666666667, 0, -7.09215686274510,...
     0, 54.9711779448622, 0, -529.124242424242,...
     0, 6192.12318840580, 0, -86580.2531135531,...
     0, 1425517.16666667, 0, -27298231.0678161,...
     0, 601580873.900642, 0, -15116315767.0922,...
     0, 429614643061.167, 0, -13711655205088.3,...
     0, 4.88332318973593e+14, 0, -19296579341940068];
  bn = length(b);


  asx = 10;			% value of x beyond which asymptotic is
				% believed to work

  if (~isreal(x)) error('x must be real.'); end;
  if (length(n)~=1) usage('n must be a scalar.'); end;
  n=floor(n);
  if (n==0)			% calling psi(x) for digamma
    f=psi(x);
    return;
  end;



  [row, col] = size(x);		% stroing dimensions of x
  x=x(:);			% restating x as a vector

  
  xx=max(x, asx+ (x-floor(x)));	% the value to go into the asymptotic formula
  recur= max(ceil(asx-x), 0);	% number of recursions needed to get to
				% that value

  oxx  = 1./(xx-1);		% temporary variables; note that x+1
				% enters the stirling approx.
  m1n  = (-1)^n;		
  nf   = prod([1:n]);
  nm1f = nf/n;
  f = m1n*(-1)*nm1f.*oxx.^n + 0.5*nf*m1n.*oxx.^(n+1);
  for i=2:2:bn
    % f +=  (-1)*m1n*prod([i:i+n-1])/i*b(i).*oxx.^(i+n);
    f = f + (-1)*m1n*prod([i:i+n-1])/i*b(i).*oxx.^(i+n);
  end;

  %% accounting for recursion that brought x's to asymptotic regime
  for i=1:length(x)
    if(recur(i))
      % f(i)+= m1n*(-1)*nf*sum((x(i)+ [0:1:(recur(i)-1)]).^(-n-1));
      f(i)=f(i) + m1n*(-1)*nf*sum((x(i)+ [0:1:(recur(i)-1)]).^(-n-1));
    end;
  end;
    
  f=reshape(f,row,col);		% reshaping to original shape
  


