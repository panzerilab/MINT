function B = B_xiK (xi, K)
  %% Usage:  B = B_xiK (xi, K)
  %%
  %% Calculates the value of B (K*\beta) that corresponds to the 
  %% a priori expected value of entropy equal to xi. Spline 
  %% interpolation is involved. To speed up generation, the function 
  %% stores the spline tables once generated in the files 
  %% bxi/BxisplinesXXX.mat, where XXX is the number of bins, K.
  %%
  %% Inputs:
  %%    xi -- value of the entropy, any dimensional object.
  %%    K  -- number of bins, scalar.
  %%
  %% Output:
  %%    B  -- value of B that corresponds to xi; same dimensionality as
  %%          xi.
  %%
  %% Depends on:
  %%    xi_KB.m
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

  
  maxcounter = 200;  % 40;

  if (length(K) ~=1) 
    error('Wrong dimension of K.');
  end;

  if (any(xi >log(K)))
    error('Too large xi -- bigger than log(K).');
  end;
      
  if (any(xi <0))
    error('Too small xi -- smaller than 0.');
  end;
  
  %% define the interp variables and see if they exist
  %% in memory already for the correct value of K
  global K_interp_in
  global Bxi_interp_in
  

  %% if no K_interp_in defined, loaded =0
  eval(['loaded = length(K_interp_in);'],'loaded=0;');
  %% if loaded, check if the correct one is loaded
  if(loaded)
    loaded = (K_interp_in==K);
  end;
  
  %% we will have an interpolation table for B=B(xi). 
  %% we will lookup values in this table, then polish them
  %% with Newton-Raphson. It is important, that lookup gets the 
  %% initial guess for B which is smaller than actual B, because 
  %% dxi/dB is small at large B and we can easily overshoot
  %% to negative B if the original B was too large.
 if(~loaded)			% no variable in memory
    disp(['      B_XIK: Creating the B(xi) interpolation data for K=' num2str(K) '.']);
    K_interp_in = K;            % need to create it

    %% first guess at intervals for interpolation table
    step=1e-2;
    b1  = 10;			% first crossover
    b2  = 100*K;		% second crossover
    b=[0:step:b1 b1*exp([step:step:log(b2)])]; 
%(10+step*K:step*K:10*K K*(10+step*K)*exp(0:2*step:20)];
    %% create the interpolation table
    Bxi_interp_in = [xi_KB(K,b);b];
    %% now using the crude interpolation table put uniform spacings
    %% onto xi and find corresponding values for B, and put both
    %% into the new interpolation table
    dxi = [1e-3:log(K)*1e-3:log(K)-1e-3];
    Bxi_interp_in = [dxi; B_xiK(dxi,K)];
    disp('      B_XIK: Done creating interpolation data.');
 end;


  %% now the interpolation table is in the memory, get 
  %% the (approximate) values of B given xi
  [row, col] =size(xi);
  xi=xi(:);
  B=Bxi_interp_in(2,lookup(Bxi_interp_in(1,2:size(Bxi_interp_in,2)-1), xi)+1);
  B=B(:);  
  
  
  %% now do the Newton-Raphson polishing to find B
  %% the absolut error of 1e-10 is needed for the iterations to converge
  %% near xi=0; these values will later be repplaced by the asymptotic
  %% expansion
  counter = 0;

  F   =  9999999;
  %% again, checking precision in xi, not in B
  while ~(all(abs(F) < abs(xi*1e-13 +1e-13)) | (counter>maxcounter));
     counter = counter + 1;
     F   =  (xi_KB(K,B))-xi;
     dF  =  dxi_KB(K,B);
     dB  =  - F./dF;
     dB(abs(F)<(xi*1e-13)) = 0;	% killing the problem with the precision
                                % loss; effectively we converge to 1e-12 
				% in xi , not B
    B  = B + dB;
  end; 
  if(counter>=maxcounter)
    error(['Newton-Raphson root polishing did not converge after ' ...
	     int2str(counter) ' iterations. Problems are likely.']);
  end;
  
  %% make replacement of the asymptotically large beta
  repl_large = find((B>100*K) & (xi~=log(K)));
  B(repl_large)= (K-1)./(2*(log(K)- xi(repl_large)));
  B(find(xi==log(K))) = Inf;

  %% make replacement of the asymptotically small beta
  repl_small = find(xi< (pi^2 * 1e-6 /6));
  B(repl_small)= (1+1/(K-1))*xi(repl_small)*(6/pi^2);
  

  B=reshape(B,row,col);  
  
%%%%%%%%%%%%%%%%%
%%% old function version


% function B = B_xiK (xi, K)
% 
%   %% Usage:  B = B_xiK (xi, K)
%   %%
%   %% Calculates the value of B (K*\beta) that corresponds to the 
%   %% a priori expected value of entropy equal to xi. Spline 
%   %% interpolation is involved. To speed up generation, the function 
%   %% stores the spline tables once generated in the files 
%   %% bxi/BxisplinesXXX.mat, where XXX is the number of bins, K.
%   %%
%   %% Inputs:
%   %%    xi -- value of the entropy, any dimensional object.
%   %%    K  -- number of bins, scalar.
%   %%
%   %% Output:
%   %%    B  -- value of B that corresponds to xi; same dimensionality as
%   %%          xi.
%   %%
%   %% Depends on:
%   %%    xi_KB.m
%   %%  
%   %%  (c) Ilya Nemenman, March 2002
%   
%   maxcounter = 40;
% 
%   if (length(K) ~=1) 
%     error('Wrong dimension of K.');
%   end;
% 
%   if (any(xi >log(K)))
%     error('Too large xi -- bigger than log(K).');
%   end;
%       
%   if (any(xi <0))
%     error('Too small xi -- smaller than 0.');
%   end;
%   
%   %% define the interp variables and see if they exists
%   %% in memory already for the correct value of K
%   global K_interp_in
%   global Bxi_interp_in
% 
%   %% if no K_interp_in defined, loaded =0
%   eval(['loaded = length(K_interp_in);'],'loaded=0; ');
%   %% if loaded, check if the correct one is loaded
%   if(loaded)
%     loaded = (K_interp_in==K);
%   end;
%   
%   %% we will have an interpolation table for B=B(xi). 
%   %% we will lookup values in this table, then polish them
%   %% with Newton-Raphson. It is important, that lookup gets the 
%   %% initial guess for B which is smaller than actual B, because 
%   %% dxi/dB is small at large B and we can easily overshoot
%   %% to negative B if the original B was too large.
%   if(~loaded)			% no variable in memory
%     disp(['      B_XIK: Creating the B(xi) interpolation data for K=' num2str(K) '.']);
%     K_interp_in = K;            % need to create it
% 
%     %% first guess at intervals for interpolation table
%     step=1e-2;
%     b1  = 10;			% first crossover
%     b2  = 100*K;		% second crossover 
%     b=[0:step:b1 b1*exp([step:step:log(b2)])] 
% %(10+step*K:step*K:10*K K*(10+step*K)*exp(0:2*step:20)];
%     %% create the interpolation table
%     Bxi_interp_in = [xi_KB(K,b);b]
%     %% now using the crude interpolation table put uniform spacings
%     %% onto xi and find corresponding values for B, and put both
%     %% into the new interpolation table
%     dxi = [1e-3:log(K)*1e-3:log(K)-1e-3];
%     Bxi_interp_in = [dxi; B_xiK(dxi,K)];
%     disp('      B_XIK: Done creating interpolation data.');
%   end; % 
% 
% 
%   %% now the interpolation table is in the memory, get 
%   %% the (approximate) values of B given xi
%   [row, col] =size(xi);
%   xi=xi(:);
%   %
%   Bxi_interp_in(1,2:size(Bxi_interp_in,2))
%   % table = Bxi_interp_in(1,2:size(Bxi_interp_in,2)-1)
%   B=Bxi_interp_in(2, lookup(Bxi_interp_in(1,2:size(Bxi_interp_in,2)-1), xi)+1 );
%  % B=Bxi_interp_in(2,lookup(Bxi_interp_in(1,2:size(Bxi_interp_in,2)-1), xi) +1);
%   B=B(:);
% 
% 
%   %% now do the Newton-Raphson polishing to find B
%   %% the absolut error of 1e-10 is needed for the iterations to converge
%   %% near xi=0; these values will later be repplaced by the asymptotic
%   %% expansion
%   counter = 0;
%   for counter=counter+1;
%       F   =  xi_KB(K,B)-xi;
%       dF  =  dxi_KB(K,B);
%       dB  =  - F./dF;
%       while ~((all(abs(F) < abs(xi*1e-13 +1e-13)) | (counter > maxcounter)));
%           % do
%           % counter++;
%           counter = counter + 1;
%           F   =  xi_KB(K,B)-xi;
%           dF  =  dxi_KB(K,B);
%           dB  =  - F./dF;
%           %    if(length(xi)==1) [xi F B dB], endif
%           dB(abs(F)<(xi*1e-13)) = 0;	% killing the problem with the precision
% 				% loss; effectively we converge t 1e-12
% 				% in xi , not B
%                 B  = B + dB;
%     
%        end;
%    end;
%    % B  += dB;
%     %% again, checking precision in xi, not in B
%   % until (all(abs(F) < abs(xi*1e-13 +1e-13)) | (counter>maxcounter));
%   
% 
%   if(counter>=maxcounter)
%     error(['Newton-Raphson root polishing did not converge after ' ... %\
% 	     int2str(counter) ' iterations. Problems are likely.']);
%   end;
%       
% %     end;
% %  end;
%   
%   
%   %% make replacement of the asymptotically large beta
%   repl_large = find((B>100*K) & (xi~=log(K)));
%   B(repl_large)= (K-1)./(2*(log(K)- xi(repl_large)));
%   B(find(xi==log(K))) = Inf;
% 
%   %% make replacement of the asymptotically small beta
%   repl_small = find(xi< (pi^2 * 1e-6 /6));
%   B(repl_small)= (1+1/(K-1))*xi(repl_small)*(6/pi^2);
%   
% 
%   B=reshape(B,row,col);  
% 
% % endfunction    
% 
