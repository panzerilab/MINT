function [S_as dS_as D] = asymp_val (path, inname, ...)
  %% Usage:  [S_as, dS_as, D] = 
  %%               asymp_val (path, inname, ... )
  %%
  %% The function parses a given set of preprocessed data files 
  %% with occurence statistics of words of different length and
  %% calculates entropies of distribtions of these words using 
  %% asymptotic formulas in Nemenman, 2002.
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
  %% Optional arguments:
  %%     fstart - scalar integer, start processing with this word length;
  %%              default value -- 1;
  %%     fend   - scalar integer, stop processing at this word length.
  %%              default value -- process while data is available.
  %%
  %% Out:
  %%     Vectors of calculated entropy estimates are both returned, and 
  %%     also stored in the file with name 
  %%        [inname + 'entr_as.txt']
  %%     in the directory 'path'. The length of all return variables is 
  %%     determined by how many different data files are present, and 
  %%     results for i'th data file are in i'th element of each return
  %%     variable. If mistake occured while processing some file, the
  %%     routine tries to save information calculated for smaller word 
  %%     sizes. 
  %%     S_as    - asymptotic entropy estimate;
  %%     dS_as   - the standard deviation of the estimate;
  %%     D       - number of coincidences;
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


  S_as   = zeros(1,1000);
  dS_as  = zeros(1,1000);
  D      = S_as;

  %% which files to process?
  if (nargin==2)		% process all files
    i    = 1;
    iend = Inf;
  elseif (nargin==3)		% process files from nargin(3) and up
    va_start();
    i    = floor(va_arg());
    iend = Inf;
  else				% nargin==4 or more
    va_start();			% do from nargin(3) to nargin(4)
    i     = floor(va_arg());
    iend  = floor(va_arg());
  end; 
    

  while ((i<=iend) && ...
	 ( exist(fullname =[path '/' inname '_' num2str(i) '.txt'])==2))

      disp(['ASYMP_VAL: Processing file ' fullname '.' ]);
      eval(['load ' fullname]);
      varname = [inname '_' num2str(i)];
      
      nx = eval([varname ' (1,:);']);
      kx = eval([varname ' (2,:);']);
      
      nz = find(nx >0);
      nx = nx(nz);
      kx = kx(nz);
      N  = sum(kx.*nx);
      K1 = sum(kx);
      D(i)= N-K1;
      S_as(i) = -psi(1) -log(2) +2*log(N)-psi(D(i));
      dS_as(i)= sqrt(polygamma(1,D(i)));

      
      i=i+1;

  end; 


  S_as   = S_as(1:(i-1));
  dS_as  = dS_as(1:(i-1));
  D = D(1:(i-1));

  %% saving all to a file
  outname = [path '/' inname  '_entr_as.txt'];
  eval(['save -ascii ' outname ' S_as dS_as D']);

