function [S_nsb, dS_nsb, S_cl, dS_cl, xi_cl, S_ml, errcode] = \
      entropyN (path, inname, asize, precision, qfun, ...)
  %% Usage:  [S_nsb, dS_nsb, S_cl, dS_cl, xi_cl, s_ml, errcode] = 
  %%               entropyN (path, inname, asize, precision, qfun, ...)
  %%
  %% The function parses a given set of preprocessed data files 
  %% with occurence statistics of words of different length and
  %% calculates entropies of distribtions of these words. 
  %%
  %% In:
  %%     path   - path to the data files; string.
  %%     inname - generic name of the data files, string. The name of a 
  %%              particular data file that has statistics of 
  %%              words of length i is 
  %%              [inname + '_' + i '.txt'].
  %%              Files are assumed to be created by the software 
  %%              from the module Entropy_Parsing. Their format is:
  %%              [column of nx  column of kx], where nx is the number
  %%              of counts in a bin, and kx is the number of such bins.
  %%              Example:
  %%                 0   4.08846e+61
  %%                 1   1.20610e+06
  %%                 2   50
  %%                 3   55
  %%              The files may or may not contain the enry for nx=0.
  %%     asize  - scalar, alphabet size used to create the files;
  %%     precision - scalar, requested relative tolerance for all 
  %%              calculations;
  %%     qfun   - which integration routine to use:
  %%                 1 - quad;
  %%                 2 - quadg;
  %% Optional arguments:
  %%     fstart - scalar integer, start processing with this word length;
  %%              default value -- 1;
  %%     fend   - scalar integer, stop processing at this word length.
  %%              default value -- process while data is available.
  %%
  %% Out:
  %%     Vectors of calculated entropy estimates are both returned, and 
  %%     also stored in the file with name 
  %%        [inname + 'entr.txt']
  %%     in the directory 'path'. The length of all return variables is 
  %%     determined by how many different data files are present, and 
  %%     results for i'th data file are in i'th element of each return
  %%     variable. If mistake occured while processing some file, the
  %%     routine tries to save information calculated for smaller word 
  %%     sizes. 
  %%     S_nsb   - entropy estimate by the NSB method, scalar;
  %%     dS_nsb  - the standard deviation of the estimate;
  %%     S_cl    - entropy at the saddle point;
  %%     dS_cl   - standard deviation at the saddle point;
  %%     xi_cl   - saddle point;
  %%     errcode - error code; this is build as the error code from
  %%               finding the saddle, plus 10 times the error
  %%               code of the normalization integral, 100 times
  %%               the error code of the S integral, and 1000 times
  %%               the error code of the S^2 integral. The saddle
  %%               finding error code are (see max_evidence.m) 
  %%                   0 - all ok;
  %%                   1 - no coincidences; saddle point evaluation 
  %%                       invalid (wide variance);
  %%                   2 - all data coincides; saddle point evaluation
  %%                       invalid - Bcl close to zero;
  %%                   3 - no convergence in Newton-Raphson root
  %%                       finding of B_cl.
  %%               The integration erros are (see dqag.f documentation
  %%               in QUADPACK):
  %%                   0 - all ok (but check estimated error anyway);
  %%                   1 - maximum allowed number of Gauss-Kronrod 
  %%                       subdivisions has been achieved;
  %%                   2 - roundoff error is detected; the requested 
  %%                       tolerance cannot be achieved;
  %%                   3 - extremely bad integrand behavior encountered;
  %%        	       6 - invalid input.
  %%
  %% Depends on:
  %%     max_evidence.m, integrand_1.m, integrand_S.m, integrand_S2.m
  %%     mlog_evidence.m
  %%
  %% This is the only MatLab incompatible file; it was use specifically
  %% for the fly data analysis, and it is probably not portable (and not
  %% needed) for other applications of NSB. The users are welcome to
  %% port it if needed.  
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

  S_nsb   = zeros(1,1000);
  dS_nsb  = S_nsb;
  S_cl    = S_nsb;
  dS_cl   = S_nsb;
  xi_cl   = S_nsb;
  errcode = S_nsb;
  S_ml    = S_nsb;

  %% which files to process?
  if (nargin==5)		% process all files
    i    = 1;
    iend = Inf;
  elseif (nargin==6)		% process files from nargin(5) and up
    va_start();
    i    = floor(va_arg());
    iend = Inf;
  else				% nargin==6 or more
    va_start();			% do from nargin(5) to nargin(6)
    i     = floor(va_arg());
    iend  = floor(va_arg());
  endif
    
  while ((i<=iend) && \
	 ( exist(fullname =[path '/' inname '_' num2str(i) '.txt'])==2))
    try
      disp(['ENTROPYN: Processing file ' fullname '.' ]);
      eval(['load ' fullname]);

      varname = [inname '_' num2str(i)];

      nx = eval([varname '(:,1);']);
      nx = nx';
      kx = eval([varname '(:,2);']);
      kx = kx';

      nz = find(nx >0);
      nx = nx(nz);
      kx = kx(nz);
      K  = asize^i;

      [S_nsb(i), dS_nsb(i), S_cl(i), dS_cl(i), xi_cl(i), errcode(i)] = \
	  find_nsb_entropy(kx,nx,K,precision, qfun);
      N=sum(kx.*nx);
      S_ml(i)= -sum(nx/N.*log(nx/N).*kx);

      i++;
    catch
      if (length(__error_text__))
	disp(['error: ENTROPYN: Error occured when working with '\
	      'word length of ' num2str(i) '. ' __error_text__ ]);
      else
	disp(['error: ENTROPYN: An unknown error occured when '\
	      'working with word length of ' num2str(i) '. ']);
      endif
      %% finish the cycle
      break;
    end_try_catch
  endwhile


  S_nsb   = S_nsb(1:(i-1));
  dS_nsb  = dS_nsb(1:(i-1));
  S_cl    = S_cl(1:(i-1));
  dS_cl   = dS_cl(1:(i-1));
  xi_cl   = xi_cl(1:(i-1));
  errcode = errcode(1:(i-1)); 
  S_ml    = S_ml(1:(i-1));

  %% saving all to a file
  outname = [path '/' inname  '_entr.txt'];
  eval(['save -ascii ' outname ' S_nsb dS_nsb S_cl dS_cl xi_cl S_ml errcode']);

endfunction
