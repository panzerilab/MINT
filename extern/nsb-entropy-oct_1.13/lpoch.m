function f = lpoch (x, n)
  %% Usage:  f = lpoch (x, n)
  %%
  %% Calculates the logarithm of Pochhammer's symbol:
  %% lpoch(x,n) = log (gamma(x+n)) - log(gamma(x)).
  %% For relatively small (<1e15) x lgamma.m is called.
  %% For large x a crude asymptotic expansion is calcualted.
  %%
  %% Input:
  %%   x - argument, any matrix;
  %%   n - order of the symbol, a scalar.
  %%
  %% Output:
  %%   f - value of the function, same dimensions as x
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



  if (length(n)~=1)
    usage('n must be a scalar.');
  end;
  
  
  %% very crude approximation!!!!!!!!!!!!!!!!!!
  
  [row, col] = size(x);		% storing dimensions of x
  x=x(:);			% restating x as a vector
  f = zeros(size(x));
  
  %% for asymptotic expansion we need n<x, and huge x
  asymp = find((n<(1e-8*x)) & (x>1e10));
  nasymp = find((n>=(1e-8*x)) | (x<=1e10));
  
  f(nasymp)=gammaln(x(nasymp)+n)-gammaln(x(nasymp));
  %% first term of asymptoti expansion
  f(asymp) = n*log(x(asymp)+n);
  i=0;
  if(length(asymp)>=1)
    %% now we are going to add up higher order terms which come
    %% from x*log(1+n/x). They are of the form n^(i+1)/x^i. Terms
    %% that come from higher orders in Stirling expansion are
    %% n^i/x^i, and are small; we do not take them into the account.
    delta = n;
    for i=i+1;
      while ~(all(abs(delta) < 1e-8*f(asymp)));
        i=i+1;
        delta = -delta.*(n./x(asymp)) *i /(i+1);
        max(abs(delta))
        f(asymp) = f(asymp) + delta;
        % until (all(abs(delta) < 1e-8*f(asymp)));
      end;
    end;
  end;
  
  f=reshape(f,row,col);		% reshaping to original shape



