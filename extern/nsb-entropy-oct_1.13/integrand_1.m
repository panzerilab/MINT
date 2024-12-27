function f = integrand_1 (xi)
  % Usage:  f = integrand_1 (xi)
  %
  % Computes the value of the integrand for the evidence integral 
  % including the prior over xi (see NSB method for 
  % calculating entropies of discrete pdf's).
  %
  % In:
  %     xi  - expected value of entropy, integration variable.
  % Globals:
  %     nsb_kx_quad, nsb_nx_quad - row vectors; exactly kx(i) 
  %           bins had nx(i) counts;
  %           IMPORTANT: bins with 0 counts are not indexed
  %           (that is, there is no nx(i) == 0) 
  %     nsb_K_quad - number of bins; it can be Inf if the number 
  %           of bins is unknown and assumed infinite; 
  %     nsb_mlog_quad -value of the negative log-evidence (action) 
  %           at the saddle point (will be factored from the integrand 
  %           to avoid overflows)
  % Out:
  %     f   - value of the integrand, same size as xi
  %
  % Depends on:
  %     mlog_evidence.m, prior_xi.m, B_xiK.m
  %
  % (c) Ilya Nemenman, 2002--2011
  % Distributed under GPL, version 2
  % Copyright (2006). The Regents of the University of California.

  % This material was produced under U.S. Government contract
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

  % octave does not allow to pass parameters to the 
  % integrand functions; the parameters thus have to be passed
  % to the functions as global variables
  
  global nsb_kx_quad nsb_nx_quad nsb_K_quad nsb_mlog_quad 


  B=B_xiK(xi, nsb_K_quad);

  
  f =  exp( (-mlog_evidence(B,nsb_kx_quad,nsb_nx_quad,nsb_K_quad))'  ...
  	 + (nsb_mlog_quad*ones(size(xi))) ).*(prior_xi(xi,nsb_K_quad));
  
 
