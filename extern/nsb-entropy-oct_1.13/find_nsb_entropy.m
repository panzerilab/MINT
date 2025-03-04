function [S_nsb, dS_nsb, S_cl, dS_cl, xi_cl, S_ml, errcode] = ...
    find_nsb_entropy (kx, nx, K, precision, qfun)
  % Usage:  [S_nsb, dS_nsb, S_cl, dS_cl, xi_cl, S_ml,errcode] =
  %                         find_nsb_entropy (kx, nx, K, precision,
  %                                           qfun)
  %
  %
  % The function calculates the estimate of the entropy and its
  % standard deviation for discrete probability distributions
  % by the NSB method.
  %
  % In:
  %     kx, nx - row vectors; exactly kx(i) bins had nx(i) counts;
  %              IMPORTANT: bins with 0 counts are not indexed
  %              (that is, there is no nx(i) == 0)
  %     K - number of bins; it can be Inf if the number of bins is
  %         unknown and assumed infinite;
  %     precision - relative precision for calculations;
  %     qfun    - which integration routine to use:
  %                  1 - quad;
  %                  2 - gaussq -- gaussian quadrature;
  %
  % Out:
  %     S_nsb   - entropy estimate by the NSB method, scalar;
  %     dS_nsb  - the standard deviation of the estimate;
  %     S_cl    - entropy at the saddle point;
  %     dS_cl   - standard deviation at the saddle point;
  %     xi_cl   - saddle point;
  %     S_ml    - maximum likelihood (naive) entropy estimate;
  %     errcode - error code; this is build as the error code from
  %               finding the saddle, plus 10 times the error
  %               code of the normalization integral, 100 times
  %               the error code of the S integral, and 1000 times
  %               the error code of the S^2 integral. The saddle
  %               finding error code are (see max_evidence.m)
  %                   0 - all ok;
  %                   1 - no coincidences; saddle point evaluation
  %                       invalid (wide variance);
  %                   2 - all data coincides; saddle point evaluation
  %                       invalid - Bcl close to zero;
  %                   3 - no convergence in Newton-Raphson root
  %                       finding of B_cl.
  %                   4 - other errors, possibly serious
  %
  %               If quad integrator is used, then the integration
  %               erros are (see dqag.f documentation function
  %                in QUADPACK) (this is only true for OCTAVE, not MATLAB)
  %                   0 - all ok (but check estimated error anyway);
  %                   1 - maximum allowed number of Gauss-Kronrod
  %                       subdivisions has been achieved;
  %                   2 - roundoff error is detected; the requested
  %                       tolerance cannot be achieved;
  %                   3 - extremely bad integrand behavior encountered;
  %         	      6 - invalid input.
  %               If quadg integrator is used, no integration error
  %               code is reported (only saddle finding is), and one
  %               has to track on-screen messages.
  %
  % Depends on:
  %     max_evidence.m, integrand_1.m, integrand_S.m, integrand_S2.m
  %     mlog_evidence.m, means.m
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


  % integrands to call
  integs=  ['integrand_1 ';'integrand_S ';'integrand_S2'];
  msgs  =   ['normalization';'S            ';'S^2          '];

  % integration error codes messages
  ecmsgs= ['Maximum allowed number of Gauss-Kronrod subdivisions has been achieved.', ... % error code 1
    'Roundoff error is detected; the requested tolerance cannot be achieved.', ... % error code 2
    'Extremely bad integrand behavior encountered.', ... % error code 3
    ' ', ...     % empty code 4
    ' ', ...     % empty code 5
    'Invalid input.'];	% error code 6

  % in the worst, case integrate over (almost) the whole range
  % of xi:
  edges  = [precision, log(K) - precision];
  xi_lim = zeros(1,2);

  % finding the position of the peak, variance around it
  % and entropy at the peak
  [B_cl, xi_cl, dS_cl, errcode] = max_evidence(kx, nx, K, precision);
  
  % it's impossible to pass parameters to integrands with octave's
  % quad routine. Thus we use global variables to store parameters.
  global nsb_K_quad nsb_kx_quad nsb_nx_quad nsb_mlog_quad nsb_val_quad
  
  nsb_K_quad = K;
  nsb_kx_quad= kx;
  nsb_nx_quad= nx;

  % is it MatLab or Octave that is being run? MatLab version() has a 
  % release (R) number in it
  MLab = length(findstr(version,'R'));
  if(~MLab)
    quad_options('abs',0);
    quad_options('rel',precision);
  end

  if((errcode==1) | (errcode==2) | (errcode==3) | (errcode==4))
    % if no coincidences, or no covergence in NR polishing,
    % or bad power-series expansion results
    % disp(['  FIND_NSB_ENTROPY: Switching to integration over the full range of xi due to ' ...
    %     'prior errors.'])
    if (B_cl>0)
        nsb_mlog_quad  = mlog_evidence(B_cl, kx, nx, K); % this is the best guess value for mlog_evidence
                    % but it can be wrong if finding the saddle point
                    % worked bad; still -- better than nothing
        S_cl = meanS(B_cl, kx, nx, K); % the best guess for S_cl but maay be wrong
    else
        nsb_mlog_quad=0; %a real bad error when Bcl is negative
        S_cl=NaN;
    end
    xi_lim(1)  = edges(1);
    xi_lim(2)  = edges(2);
    delta = 0;
    dS_cl = NaN	;	
    xi_cl = Inf ;
  else
    S_cl = meanS(B_cl, kx, nx, K);
    % disp(['  FIND_NSB_ENTROPY: Expect S to be near ' num2str(S_cl) ', and sigma(S) near ' ...
    %     num2str(dS_cl) '.']);
    nsb_mlog_quad = mlog_evidence(B_cl, kx, nx, K); % value at the saddle

    % limits of integration determination
    % delta is the interval around the peak on which
    % gaussian approximation falls to precision
    % sqrt(2) is there because erf is defined as
    % int_0^x exp(-t^2)dt
    delta = erfinv(1-precision)*sqrt(2);
    % disp(['  FIND_NSB_ENTROPY: Integrating around the peak.']);
  end; 


  nsb_val_quad  = zeros(1,size(integs,1));
  ec   = nsb_val_quad;
  nfun = nsb_val_quad;
  err  = nsb_val_quad;
  for i=1:size(integs,1);
    show = size(integs,1);
    % if integrating around the peak, need to find the limits of
    % integration
    if(delta>0)
        cent = feval(deblank(integs(i,:)),xi_cl);
        % if the integral was purely Gaussian, the value at +-delta
        % would've been exp(-(delta^2)/2), and we would've achieved
        % the required precision. Just in case the integral is
        % non-gaussian, we require the value at the limits of integration
        % to be slightly (ten times) less than that.
        good_enough = 0.1*cent*exp(-(delta^2)/2);
        % increase the integration window, until the integrand is small
        % at the edges; start with xi_cl +- delta*dxi;
        % the integration window will be (delta+wnd(1))*dS_cl on the
        % left, and (delta+wnd(2))*dS_cl on the right; we will
        % increase wnd in steps of .5
        wnd = [-1,-1];
        worst = [1, 2];		% both limits are bad
        best  = [1, 2];

	    limval(worst) = good_enough*2;
	    limval(best) = good_enough*2;
        finished=0;
        while(~finished)
            % if worst is at egdes, change the best limit
            if any(xi_lim(worst) == edges);
                % the *=1.2 clause is added so that we never get stuck with
                % the limits growing too slowly
                if (wnd(best)>10)  
                  wnd(best)=wnd(best)*1.2 ; 
                else
                  wnd(best) = wnd(best) + 0.5;
                end;
            else
                % else change the worst limit
                if(wnd(worst)>10) 
                  wnd(worst)=wnd(worst)*1.2;
                else
                  wnd(worst) = wnd(worst) + 0.5;
                end;
            end;
            xi_lim(1) = max(edges(1), xi_cl - (delta+wnd(1))*dS_cl);
            xi_lim(2) = min(edges(2), xi_cl + (delta+wnd(2))*dS_cl);
            
            % values of the integrand at the limits
            limval(worst) = feval(deblank(integs(i,:)),xi_lim(worst));
            % we don't need  the next line, except when 'worst' is an edge
            % however, it's easier to put the line in and forget
            % about this extra unnecessary evaluation taking place
            % most of the time
            limval(best) = feval(deblank(integs(i,:)),xi_lim(best));
            % the worst and the best limits
            [tmp, worst] = max(limval);
            [tmp, best]  = min(limval);
            % we end the cycle when either the worst of the limits is good
            % enough, or the best limit is good enough and the worst limit
            % is one of the edges, or both limits are at the edges
            if (limval(worst)<good_enough)
              finished=1;
            end;
            if (limval(best)<good_enough && any(xi_lim(worst)==edges))
              finished=1;
            end;
            if  all(xi_lim==edges)
              finished=1;
            end;
        end; 
    end;


    % disp(['  FIND_NSB_ENTROPY: Doing ' deblank(msgs(i,:)) ' integral. ' ...
    %     'Limits: ' num2str(xi_lim(1)) ' < xi < ' num2str(xi_lim(2)) ...
    %     '.']);

    if(qfun==1)
      % quad integrator is used
      if(MLab)      % if being executed by matlab
	  % notice that "precision" here is an absolute error, not 
	  % the relative one; however, we expect the integrand to be ~1
	  % since we did the scale extraction through mlogevidence;
	  % however, this may still remain a problem
	  nsb_val_quad(i) = quad(deblank(integs(i,:)), xi_lim(1), ...
				 xi_lim(2) , precision,0);
	  ec(i)=0;
	  err(i)=precision;
      else        % octave execution 
	      [nsb_val_quad(i), ec(i), nfun(i), err(i)] = ...
              quad (deblank(integs(i,:)), xi_lim(1), xi_lim(2), ... 
		      [0  precision]);
      end
    elseif(qfun==2)
      % gaussq integrator is used, no error tracking here; same calls
      % for MatLab and Octave
      [nsb_val_quad(i) tol] = gaussq (deblank(integs(i,:)), xi_lim(1), xi_lim(2), precision,[]);
      ec(i) = 0;		% no error reporting
      err(i) = tol*nsb_val_quad(i); %estimated relative error
    else
      % wrong qfun
      error('Invalid choice for the integrator routine.');
    end;

    % now do the error tracking
    if(ec(i))
      % disp(['warning: FIND_NSB_ENTROPY: Problem in ' deblank(msgs(i,:)) ...
      %       ' integral. ' deblank(ecmsgs(ec(i),:))]);
    end
    if(err(i)>precision)
      % disp(['warning: FIND_NSB_ENTROPY: Precision of ' num2str(precision) ...
      %       ' required, but only ' num2str(err(i)) ' achieved.']);
      % disp(['                           ' num2str(nfun(i)) ' function ' ...
      %       'evaluations perfomed.']);
    end
  end


  S_nsb  = nsb_val_quad(2)/nsb_val_quad(1);
  dS_nsb = sqrt(abs(nsb_val_quad(3)/nsb_val_quad(1)- S_nsb^2));
  % disp(['  FIND_NSB_ENTROPY: Found S = ' num2str(S_nsb) ', and sigma(S) = ' ...
	% num2str(dS_nsb) '.'])
  errcode =  errcode + ec(1)*10 + ec(2)*100 + ec(3)*1000;
  
  N  = sum (kx.*nx);
  S_ml = -sum(nx/N.*log(nx/N).*kx);
  
  clear nsb_K_quad nsb_kx_quad nsb_nx_quad 


  
  
