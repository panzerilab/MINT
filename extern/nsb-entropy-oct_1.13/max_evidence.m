function [Bcl, xicl, dxi, errcode] = max_evidence(kx, nx, K, precision)
  % Usage:  [Bcl, xicl, dxi, errcode] = max_evidence (kx, nx, K, precision)
  %
  % The function finds the position of the minimum of the a posteriori
  % evidence and the variance around it. The integration variable is
  % \xi -- the a-priori value of the entropy.
  %
  % In:
  %     kx, nx - row vectors; exactly kx(i) bins had nx(i) counts;
  %              IMPORTANT: bins with 0 counts are not indexed
  %              (that is, there is no nx(i) == 0)
  %     K - number of bins; it can be Inf if the number of bins is
  %         unknown and assumed infinite;
  %     precision - relative precision for calculations;
  % Out:
  %     Bcl     - the classical value of B;
  %     xicl    - the classical value of xi;
  %     dxi     - std. dev. near the classical value;
  %     errcode - error code; 0 - all ok;
  %                           1 - no coincidences; saddle point
  %                               evaluation invalid (wide variance);
  %                           2 - all data coincides; saddle point
  %                               evaluation invalid - Bcl close to zero;
  %                           3 - no convergence in Newton-Raphson root
  %                               finding of Bcl;
  %                           4 - Negative increment in Newton-Raphson
  %                               polishing; the procedure is unstable
  %                           5 - error in Newton-Raphson adjustment
  %                               for fractional countsendif
  %
  % Depends on:
  %     psi.m, polygamma.m
  %
  % (c) Ilya Nemenman, 2002--2006
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


  % Ki (i=1,2...) -- number of bins that had at least i counts in them
  
  maxcounter = 200;		% maximum number of Newton iterations
  errcode = 0;			% will reassign later if mistake happens

  % we require that none nx(i)=0 so that there's no problems with
  % overflow in the case of too many bins with zero count
  if any(nx<=0) error('Bins with zero count encountered'); end;

  N  = sum (kx.*nx);		% total number of samples
  if (N<=1) error(['Too few data samples: N=' int2str(N)]); end;

  ng1= find(nx>1);		% occupancy which is greater than one
  K1 = sum (kx);		% number of bins with nozero occupancy
  K2 = sum (kx(ng1));		% number of bins with coincidences

  % for fractional occupancies
  if (all(round(nx)==nx))	% all integer occupancies
    ints = 1;
  else				% fractional occupancies
    ints = 0;
    K1i  = K1;			% renormalizing K1
    K1   = sum(kx((nx>N/K) | (nx>=1)));
    n01 = find((nx<=1)&(nx>0));	% marginal (sparsely populated) bins
  end
  ep = (N-K1)/N;		% epsilon

  if (K1>=N)			% no coincidences
    Bcl = Inf;
    xicl = Inf;
    errcode = 1;
    dxi = NaN;
    disp('warning: MAX_EVIDENCE: No coincidences.');
  elseif (K1 <= 1)		% all data coincides
    Bcl = 0;
    xicl= 0;
    errcode = 2;
    dxi = NaN;
    disp('warning: MAX_EVIDENCE: All data coincide.');
  else				    % some non-trivial value of B0, Bcl has to be
                        % calcuated
    order = 10;			% we will calculate B to this order in ep
    b=zeros(1,12);		% coefficients of the expansion of B in
                        % powers of \epsilon

    % the power series coefficients calculated by mathematica
    b(1) = (-1 + N)/(2*N);
    b(2) = (-2 + N^(-1))/3;
    b(3) = (2 + N - N^2)/(9*N - 9*N^2);
    b(4) = (2*(2 - 3*N - 3*N^2 + 2*N^3))/(135*(-1 + N)^2*N);
    b(5) = (4*(22 + 13*N - 12*N^2 - 2*N^3 + N^4))/(405*(-1 + N)^3*N);
    b(6) = (4*(-40 + 58*N + 65*N^2 - 40*N^3 - 5*N^4 + ...
        2*N^5))/(1701*(-1 + N)^4*N);
    b(7) = (4*(-9496 - 6912*N + 5772*N^2 + 2251*N^3 - 1053*N^4 - 87*N^5 ...
        + 29*N^6))/(42525*(-1 + N)^5*N);
    b(8) = (16*(764 - 1030*N - 1434*N^2 + 757*N^3 + 295*N^4 - 111*N^5 - ...
        7*N^6 + 2*N^7))/(18225*(-1 + N)^6*N);
    b(9) = (16*(167000 + 142516*N - 108124*N^2 - 66284*N^3 + 26921*N^4 + ...
        7384*N^5 - 2326*N^6 - 116*N^7 + 29*N^8))/(382725*(-1 + ...
        N)^7*N);
    b(10)= (16*(-17886224 + 22513608*N + 37376676*N^2 - 17041380*N^3 - ...
        11384883*N^4 + 3698262*N^5 + 846930*N^6 - 229464*N^7 - ...
        9387*N^8 + 2086*N^9))/(37889775*(-1 + N)^8*N);
    b(11)= (16*(-4166651072 - 3997913072*N + 2783482560*N^2 + ...
        2290151964*N^3 - 803439834*N^4 - 395614251*N^5 + ...
        108055443*N^6 + 20215218*N^7 - 4805712*N^8 - 165395*N^9 ...
        + 33079*N^10))/(795685275*(-1 + N)^9*N);
    b(12)= (32*(52543486208 - 62328059360*N - 118489458160*N^2 + ...
        47185442088*N^3 + 44875379190*N^4 - 12359832987*N^5 - ...
        5400540075*N^6 + 1272974916*N^7 + 200644800*N^8 - ...
        42495955*N^9 - 1255067*N^10 + ...
        228194*N^11))/(14105329875*(-1 + N)^10*N);


    % calculating the value of B0 as a series exansion in powers of
    % epsilon (ep)
    B0ep= N*sum(ep.^[-1:order].*b);
    if (B0ep<0)
        B0ep=precision;
        % can still recover precision with NR root finding
        disp(['warning: MAX_EVIDENCE: Series expansion for B_0 '...
            'did not converge. ']);
    end;

    % using Newton-Raphson to polish the value of B0ep and find
    % B0 to the desired precision.
    % The equation we are solving is:
    %   K1/B0 - \Sum_{I=0}^{N-1} 1/(B0+I) \Equiv F(B0) = 0
    % The derivative is:
    %   dF/dB0 = - K1/B0^2 + \Sum_{I=0}^{N-1} 1/(B0+I)^2
    B0  = B0ep;
    dB0 = 99999999999;
    counter = 0;
    while ~((abs(dB0) < abs(B0*precision)) | (counter>maxcounter));
      counter = counter + 1 ;   
      F   =  K1/B0 + psi(B0) - psi(B0+N);
      dF  =  - K1/B0^2 + polygamma(1, B0) - polygamma(1,B0+N);
      dB0 =  - F/dF;
      B0 = B0 + dB0;
    end;

    if(counter>maxcounter)
      errcode = 3;
      disp(['warning: MAX_EVIDENCE: Newton-Raphson ' ...
            'search for B_0 did not converge after ' ...
            int2str(counter) ' iterations.']);
    end;

    Bcl = B0;		% take the calculated value as the first
                    % approximation for Bcl

    order_K = 4;	% number of terms in the series, orders
                    % up to 10 are in the file other_orders_K.m

    B   = zeros(1,order_K);

    % temporary variables
    EG = - psi(1);		% euler's gamma
    pg1B0  = polygamma(1, B0);
    pg1NB0 = polygamma(1, N+B0);
    denum  = K1/B0^2 - pg1B0 + pg1NB0; % denumerator
    pg2B0  = polygamma(2, B0);
    pg2NB0 = polygamma(2, N+B0);
    pg21   = polygamma(2,1);
    pg3B0  = polygamma(3, B0);
    pg3NB0 = polygamma(3, N+B0);
    pg4B0  = polygamma(4, B0);
    pg4NB0 = polygamma(4, N+B0);

    f0   = sum(kx(ng1).*psi(nx(ng1)));
    d1f0 = sum(kx(ng1).*polygamma(1, nx(ng1)));
    d2f0 = sum(kx(ng1).*polygamma(2, nx(ng1)));
    d3f0 = sum(kx(ng1).*polygamma(3, nx(ng1)));

    % calcuating the expansion orders
    B(1) = (B0^2*(EG*K2 + f0)) / (B0^2*denum);

    B(2) = (K2*pi^2*B0 - (6*K1*B(1)^2)/B0^3 - 3*B(1)^2*pg2B0 + ...
        3*B(1)^2*pg2NB0 - 6*B0*d1f0)/(-6*denum);

    B(3) = (K2*pi^2*B(1) + (6*K1*B(1)^3)/B0^4 -(12*K1*B(1)*B(2))/B0^3 + ...
        3*K2*B0^2*pg21 - 6*B(1)*B(2)*pg2B0 + 6*B(1)*B(2)*pg2NB0 - ...
        B(1)^3*pg3B0 + B(1)^3*pg3NB0 - 6*B(1)*d1f0 - 3*B0^2*d2f0)/ ...
        (-6*denum);

    B(4) = -(-(K2*pi^4*B0^3)/90 + (K1*B(1)^4)/B0^5 - (K2*pi^2*B(2))/6 - ...
        (3*K1*B(1)^2*B(2))/B0^4 + (K1*B(2)^2)/B0^3 + ...
        (2*K1*B(1)*B(3))/B0^3 - K2*B0*B(1)*pg21 + ((B(2)^2 + ...
        2*B(1)*B(3))*pg2B0)/2 ...
        - ((B(2)^2 + 2*B(1)*B(3))*pg2NB0)/2 + ...
        (B(1)^2*B(2)*pg3B0)/2 - (B(1)^2*B(2)*pg3NB0)/2 + ...
        (B(1)^4*pg4B0)/ 24 - (B(1)^4*pg4NB0)/24 +  B(2)*d1f0 + ...
        B0*B(1)*d2f0 + (B0^3*d3f0)/6)/(-denum);

    % adding up the calculated expansion terms
    if(ints==1)
        Bcl = Bcl + sum(B.* K.^(-[1:order_K]));
    end				% if noninteger data present, do not sum
                    % the series

    % at this point K1 must return to unrenormalized value
    if (ints==0) K1=K1i; end;



    % using Newton-Raphson to polish the value of B0ep and find
    % B0 to the desired precision.
    counter = 0;
    dBcl = 999999999999;
    while ~((abs(dBcl) < abs(Bcl*precision)) | (counter>maxcounter));
      counter=counter+1;
      if(ints==1)		% integer counts
        F    =  1/K*sum(kx(ng1).*psi(nx(ng1) + Bcl/K)) - K2/K*psi(1+Bcl/K) + ...
            K1/Bcl + psi(Bcl) - psi(Bcl+N);
        dF   =  1/(K^2)*sum(kx(ng1).*polygamma(1, nx(ng1) + Bcl/K)) - ...
            K2/(K^2)*polygamma(1,1+Bcl/K) - K1/Bcl^2 + polygamma(1, Bcl) - ...
            polygamma(1, Bcl +N);
      else			% fractional counts
        F    =  1/K*sum(kx.*psi(nx + Bcl/K)) - K1/K*psi(Bcl/K) + ...
            psi(Bcl) - psi(Bcl+N);
        dF   =  1/(K^2)*sum(kx.*polygamma(1, nx + Bcl/K))  - ...
            K1/(K^2)*polygamma(1,Bcl/K)  + polygamma(1, Bcl) - ...
            polygamma(1, Bcl +N);
      end

      dBcl =  - F/dF;
      if(dBcl<0) dBcl=0; warning('MAX_EVIDENCE: negative (unstable) change.');errcode=4;end;
      Bcl = Bcl +  dBcl;

  end;         
                
  if(counter>maxcounter)
    errcode = 3;
    warning(['Newton-Raphson search for Bcl did not converge after ' ...
             int2str(counter) ' iterations.']);
    end;

    if ((errcode == 3) & (counter<=maxcounter))
      % earlier problem with NR integration for B0
      % was overcome in integration for maxcounter
      disp(['success: MAX_EVIDENCE: Recovered from previous errors ' ...
            'in Bcl determination.'])
      errcode = 0;
    end;

    if(ints==1)			% integer data
      dBcl = 1/K^2 *sum(kx(ng1).*polygamma(1, nx(ng1) + Bcl/K)) - ...
          K2/K^2*polygamma(1, 1+Bcl/K) - K1/Bcl^2 + polygamma(1, Bcl) - ...
          polygamma(1, Bcl +N);
    else			% fractional counts
      dBcl   =  1/(K^2)*sum(kx.*polygamma(1, nx + Bcl/K)) - ...
          K1/(K^2)*polygamma(1,Bcl/K)  + polygamma(1, Bcl) - ...
          polygamma(1, Bcl +N);
    end

    xicl= xi_KB(K,Bcl);
    dxi = 1/sqrt(-dBcl/dxi_KB(K,Bcl)^2);

  end


 



