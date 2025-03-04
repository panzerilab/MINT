function f = psi_asymp (x)
  %% Usage:  f = psi_aymp (x)
  %%
  %% Subleading asymptotic behavior of Psi (or Digamma) function valid 
  %% for positive argument.
  %%
  %%                       d
  %%        psi_asymp(x) = --log(Gamma(x)) - log(x)
  %%                       dx
  %%
  %% We aim at extremely large x, so we use the Stirling form for 
  %% psi, and not the Lancocs one. The stirling form has the log(x) 
  %% term in it, and subtraction of the logarithm can be done with no
  %% loss of precision.
  %%    psi_asymp(x) = psi(x) - log(x) = -1/2x 
  %%              - ...sum_{i=1}^{...infty} B_{2i} /(2i*x^{2i}),
  %%
  %% In:
  %%    x - value of the argument, any matrix.
  %% Out:
  %%    f - value of the function, same size as x.
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

  if (~isreal(x)) usage('x must be real.'); end;

  [row, col] = size(x);		% stroing dimensions of x
  x=x(:);			% restating x as a vector

  
  xx=max(x, asx+ (x-floor(x)));	% the value to go into the asymptotic formula
  recur= max(ceil(asx-x), 0);	% number of recursions needed to get to
				% that value

  oxx  = 1./xx;			% temporary variables; note that x 
				% (not x+1) enters the form of stirling we use
  oxxi = oxx.^bn;
  xx2  = xx.*xx;

  f = zeros(size(x));
  %% we start summation from the smallest term to avoid 
  %% precision loss
  for i=bn:-2:2
    f = f -  b(i).*oxxi/i;
   %  oxxi .*= xx2;
    oxxi = oxxi.*xx2;
  end;
  f = f - 0.5.*oxx; 

  %% accounting for recursion that brought x's to asymptotic regime
  for i=1:length(x)
    if(recur(i))
      f(i) = f(i) +  -sum(1./(x(i)+ [0:1:(recur(i)-1)]));
    end;
  end;
    
  f=reshape(f,row,col);		% reshaping to original shape
  


