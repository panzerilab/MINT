function [int, tol1,ix]= gaussq(fun,xlow,xhigh,tol,trace,p1,p2,p3,p4,p5,p6,p7,p8,p9)
%GAUSSQ Numerically evaluates a  integral using a Gauss quadrature.
%       The Quadrature integrates a (2m-1)th order  polynomial exactly 
%       and the  integral is of the form
%              b                              
%             Int (p(x)* Fun(x)) dx 
%              a                 
%
% CALL:
% [int, tol] = gaussq('Fun',xlow,xhigh,[reltol wfun],[trace,gn],p1,p2,....)
%
% [int, tol] = gaussq('Fun',xlow,xhigh,[reltol wfun],[trace,gn],alpha,p1,p2,....)
%
%  [int, tol] = gaussq('Fun',xlow,xhigh,[reltol wfun],[trace,gn],alpha,beta,p1,p2,....)
%
%       int = evaluated integral
%       tol = absolute tolerance abs(int-intold)
%       Fun = string containing the name of the function or a directly given 
%              expression enclosed in parenthesis. 
%xlow,xhigh = integration limits
%    reltol = relative tolerance default=1e-3
%      wfun = weight function
%         1  p(x)=1                       a =-1,   b = 1 Legendre (default)
%         2  p(x)=1/sqrt((x-a)*(b-x)),    a =-1,   b = 1 Chebyshev of the
%                                                             first kind
%         3  p(x)=sqrt((x-a)*(b-x)),      a =-1,   b = 1 Chebyshev of the 
%                                                            second kind
%         4  p(x)=sqrt((x-a)/(b-x)),      a = 0,   b = 1
%         5  p(x)=1/sqrt(b-x),            a = 0,   b = 1
%         6  p(x)=sqrt(b-x),              a = 0,   b = 1
%         7  p(x)=(x-a)^alpha*(b-x)^beta  a =-1,   b = 1 Jacobi 
%                                     alpha, beta >-1 (default alpha=beta=0)
%         8  p(x)=x^alpha*exp(-x)         a = 0,   b = inf generalized Laguerre
%         9  p(x)=exp(-x^2)               a =-inf, b = inf Hermite
%        10  p(x)=1                       a =-1,   b = 1 Legendre (slower than 1)
%
%   trace = for non-zero TRACE traces the function evaluations 
%              with a point plot of the integrand.
%      gn = number of base points and weight points to start the 
%            integration with (default=2)
%p1,p2,...= coefficients to be passed directly to function Fun:   
%                  G = Fun(x,p1,p2,...).
%
% Note that int is the common size of xlow, xhigh and p1,p2,....
% Example: a) integration of x.^2 from 0 to 2 and from 2 to 4 is done by:
%                        gaussq('(x.^2)',[0 2],[2 4])
%          b) integration of x^2*exp(-x) is done by 
%                        gaussq('(1)',0,inf,[1e-3 8],[],2)
%             or
%                        gaussq('(x.^2)',0,inf,[1e-3 8],[],0)
% See also qrule
%


%
%
%This function works just like QUAD or QUAD8 but uses a Gauss-Chebyshev
% quadrature integration scheme. The Quadrature integrates a 
%(2m-1)th order  polynomial exactly and the  integral is of the form
%                   b                                m
%                 Int (p(x)* fun(x)) dx  =  Sum ( wf_j* fun( bp_j ) )
%                  a                                j=1                   
% 
% modified version of quadc and quadg by Per A. Brodtkorb 30.03.99, 17.02.99 :
% -accept multiple integrationlimits, int is the common size of xlow and xhigh
% -optimized by only computing the integrals which did not converge.
%  - enabled the integration of directly given functions enclosed in 
%     parenthesis. Example: integration from 0 to 2 and from 2 to 4 for x is done by:
%                        gaussq('(x.^2)',[0 2],[2 4])


% Reference 
%   wfun 1: from grule.m in NIT toolbox, see ref [2] 
%   wfun 2-6: see ref [4]
%   wfun 7-10:  Adapted from Netlib routine gaussq.f see ref [1,3]
%
% [1]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [2] Davis and Rabinowitz (1975) 'Methods of Numerical Integration', page 365,
%     Academic Press.
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% [4] Abromowitz and Stegun (1954) 'Handbook of mathematical functions'


global ALPHA1 BETA1 ALPHA2
wfun=1;
if nargin<4| isempty(tol),
  tol=1e-3;
elseif length(tol)==2,
  wfun=tol(2);
   tol=tol(1);
end



istart=1; alpha1=0;beta1=0;
if (wfun==7)|(wfun==8),
  istart=2;
  if(nargin>=6&~isempty(p1)),
    alpha1=p1;
  end
end    
if wfun==7,
  istart=3;
  if nargin>=7&~isempty(p2),
    beta1=p2;
  end
end


FindWeigths=1;
switch wfun
  case 7,
    if isempty(ALPHA1)|isempty(BETA1),
    elseif ALPHA1==alpha1 & BETA1==beta1,
      FindWeigths=0;
    end
    ALPHA1=alpha1;BETA1=beta1;
      %remember what type of weights are saved as global constants
  case 8,
  if isempty(ALPHA2),
    elseif ALPHA2==alpha1,
      FindWeigths=0;
    end
    ALPHA2=alpha1;
    %remember what type of weights are saved as global constants
  otherwise,FindWeigths=0;
end


gn=2;
if nargin <5|isempty(trace) ,
  trace=0; 
elseif length(trace)==2,
  gn=trace(2);
  trace=trace(1);
end



if prod(size(xlow))==1,% make sure the integration limits have correct size
  xlow=xlow(ones(size(xhigh)));;
elseif prod(size(xhigh))==1,
  xhigh=xhigh(ones(size(xlow)));;
elseif any( size(xhigh)~=size(xlow) )
  error('The input must have equal size!')
end
[N M]=size(xlow);%remember the size of input


if any(fun=='('), %  & any(fun=='x'),
  exec_string=['y=',fun ';']; %the call function is already setup
else
  %setup string to call the function
  exec_string=['y=',fun,'(x'];
  num_parameters=nargin-5;
  for i=istart:num_parameters,
    xvar=['p' int2str(i)]; % variable # i
    if eval(['~ischar(' ,xvar,')  & length(',xvar,'(:)) ~=1' ]) ,
        if N*M==1,     
          eval(['[N M]=size(', xvar, ');']); 
        elseif  eval(['N*M~=prod(size(', xvar, '));']);
          error('The input must have equal size!')
        end     
      eval([xvar, '=' ,xvar ,'(:);']); %make sure it is a column 
      exec_string=[exec_string,',' xvar '(k,ones(1,gn) )']; %enable integration with multiple arguments
    else
      exec_string=[exec_string,',' xvar];
    end
  end
  exec_string=[exec_string,');'];
end



nk=N*M;% # of integrals we have to compute
k=(1:nk)';
int=zeros(nk,1);
tol1=int;


%gntxt=int2str(gn);% # of weights
wfuntxt= int2str(wfun);% weightfunction used
wtxt=[int2str(gn),'_',wfuntxt];% weightfunction used
eval(['global cb',wtxt,' cw',wtxt,';']);
if isempty(eval(['cb',wtxt]))|FindWeigths ,  
  % calculate the weights if necessary
  eval(['[cb',wtxt,',cw',wtxt,']=qrule(gn,wfun,alpha1,beta1);']);
end


%setup mapping parameters and execution strings
xlow=xlow(:);
jacob=(xhigh(:)-xlow(:))/2;


switch wfun,% this is clumsy and can written more tidy
  case {1 ,10},
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k);'];
  case 2,  
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2);'];
  case 3, 
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k).^2;'];
  case 4,  
    calcx_string=['(ones(nk,1),:)).*jacob(k,ones(1, gn ))*2+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k)*2;'];
  case 5,  
    calcx_string=['(ones(nk,1),:)).*jacob(k,ones(1, gn ))*2+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*sqrt(jacob(k)*2);'];
  case 6,  
    calcx_string=['(ones(nk,1),:)).*jacob(k,ones(1, gn ))*2+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*sqrt(jacob(k)*2).^3;'];
  case 7, 
    calcx_string=['(ones(nk,1),:)+1).*jacob(k,ones(1, gn ))+xlow(k,ones(1,gn ));'];
    int_string=['(ones(nk,1),:).*y,2).*jacob(k).^(alpha1+beta1+1);'];
  case {8,9}
     calcx_string=['(ones(nk,1),:));'];
     int_string=['(ones(nk,1),:).*y,2);'];
otherwise error('unknown option')
end



eval(['x=(cb',wtxt, calcx_string]); % calculate the x values            
eval(exec_string); % calculate function values  y=fun(x)                        
eval(['int(k)=sum(cw',wtxt, int_string]);       % do the integration
int_old=int;


if trace==1,
  x_trace=x(:);
  y_trace=y(:);
end


% Break out of the iteration loop for three reasons:
%  1) the last update is very small (compared to int  and  compared to tol)
%  2) There are more than 11 iterations. This should NEVER happen. 


converge='n';
for i=1:10,
  gn=gn*2;% double the # of weights
 % gntxt=int2str(gn);% # of weights
  wtxt=[int2str(gn),'_',wfuntxt];
  eval(['global cb',wtxt,' cw',wtxt]);
  if isempty(eval(['cb',wtxt]))|FindWeigths ,  
    % calculate the weights if necessary
    eval(['[cb',wtxt,',cw',wtxt,']=qrule(gn,wfun,alpha1,beta1);']);
  end
 
  eval(['x=(cb',wtxt, calcx_string]); % calculate the x values  
  eval(exec_string);                 % calculate function values  y=fun(x)
  eval(['int(k)=sum(cw',wtxt, int_string]);% do the integration


  if trace==1,
    x_trace=[x_trace;x(:)];
    y_trace=[y_trace;y(:)];
  end


  tol1(k)=abs(int_old(k)-int(k)); %absolute tolerance
  k=find(tol1 > abs(tol*int)); %| tol1 > abs(tol));%indices to integrals which did not converge
   
  if any(k),% compute integrals again
      nk=length(k);%# of integrals we have to compute again
  else
    converge='y';
    break;
  end
  int_old=int;
end


int=reshape(int,N,M); % make sure int is the same size as the integration  limits
if nargout>1,
        tol1=reshape(tol1,N,M);
end


if converge=='n',
  disp('Integral did not converge--singularity likely')
end


if trace==1,
  plot(x_trace,y_trace,'+')
end



function [bp,wf]=qrule(n,wfun,alpha,beta)
%QRULE computes abscissas and weight  factors for some selected 
%   Gaussian quadratures.   
%
%CALL:  [bp,wf]=qrule(n,wfun,alpha,beta)
%  
%  bp = base points
%  wf = weight factors
%  n  = number of base points (abscissas) (integrates a (2n-1)th order
%       polynomial exactly)
%wfun = weight function%     
%     1  p(x)=1                       a =-1,   b = 1 Legendre (default)
%     2  p(x)=1/sqrt((x-a)*(b-x)),    a =-1,   b = 1 Chebyshev of the
%                                                             first kind
%     3  p(x)=sqrt((x-a)*(b-x)),      a =-1,   b = 1 Chebyshev of the 
%                                                            second kind
%     4  p(x)=sqrt((x-a)/(b-x)),      a = 0,   b = 1
%     5  p(x)=1/sqrt(b-x),            a = 0,   b = 1
%     6  p(x)=sqrt(b-x),              a = 0,   b = 1
%     7  p(x)=(x-a)^alpha*(b-x)^beta  a =-1,   b = 1 Jacobi 
%                                     alpha, beta >-1 (default alpha=beta=0)
%     8  p(x)=x^alpha*exp(-x)         a = 0,   b = inf generalized Laguerre
%     9  p(x)=exp(-x^2)               a =-inf, b = inf Hermite
%    10  p(x)=1                       a =-1,   b = 1 Legendre (slower than 1)
%
%  The Gaussian Quadrature integrates a (2n-1)th order
%  polynomial exactly and the integral is of the form
%           b                         n
%          Int ( p(x)* F(x) ) dx  =  Sum ( wf_j* F( bp_j ) )
%           a                        j=1                          
%  See also: gaussq


% Reference 
%   wfun 1: copied from grule.m in NIT toolbox, see ref [2] 
%   wfun 2-6: see ref [4]
%   wfun 7-10:  Adapted from Netlib routine gaussq.f see ref [1,3]
%
% [1]  Golub, G. H. and Welsch, J. H. (1969)
% 'Calculation of Gaussian Quadrature Rules'
%  Mathematics of Computation, vol 23,page 221-230,
%
% [2] Davis and Rabinowitz (1975) 'Methods of Numerical Integration', page 365,
%     Academic Press.
%
% [3]. Stroud and Secrest (1966), 'gaussian quadrature formulas', 
%      prentice-hall, Englewood cliffs, n.j.
% 
% [4] Abromowitz and Stegun (1954) ''


%  By Bryce Gardner, Purdue University, Spring 1993.
% Modified by Per A. Brodtkorb 19.02.99 pab@marin.ntnu.no
% to compute other quadratures  than the default
if nargin<4|isempty(beta),
 beta=0; 
end


if nargin<3|isempty(alpha),
  alpha=0; 
end
if alpha<=-1 | beta <=-1,
  error('alpha and beta must be greater than -1')
end


if nargin<2|isempty(wfun),
  wfun=1; 
end     



switch wfun, %
  case 1,
    %  This routine computes Gauss base points and weight factors
    %  using the algorithm given by Davis and Rabinowitz in 'Methods
    %  of Numerical Integration', page 365, Academic Press, 1975.
    bp=zeros(n,1); wf=bp; iter=2; m=fix((n+1)/2); e1=n*(n+1);
    mm=4*m-1; t=(pi/(4*n+2))*(3:4:mm); nn=(1-(1-1/n)/(8*n*n));
    xo=nn*cos(t);
    for j=1:iter
      pkm1=1; pk=xo;
      for k=2:n
        t1=xo.*pk; pkp1=t1-pkm1-(t1-pkm1)/k+t1;
        pkm1=pk; pk=pkp1;
      end
      den=1.-xo.*xo; d1=n*(pkm1-xo.*pk); dpn=d1./den;
      d2pn=(2.*xo.*dpn-e1.*pk)./den;
      d3pn=(4*xo.*d2pn+(2-e1).*dpn)./den;
      d4pn=(6*xo.*d3pn+(6-e1).*d2pn)./den;
      u=pk./dpn; v=d2pn./dpn;
      h=-u.*(1+(.5*u).*(v+u.*(v.*v-u.*d3pn./(3*dpn))));
      p=pk+h.*(dpn+(.5*h).*(d2pn+(h/3).*(d3pn+.25*h.*d4pn)));
      dp=dpn+h.*(d2pn+(.5*h).*(d3pn+h.*d4pn/3));
      h=h-p./dp; xo=xo+h;
    end
    bp=-xo-h;
    fx=d1-h.*e1.*(pk+(h/2).*(dpn+(h/3).*(...
        d2pn+(h/4).*(d3pn+(.2*h).*d4pn))));
    wf=2*(1-bp.^2)./(fx.*fx);
    if (m+m) > n, bp(m)=0; end
    if ~((m+m) == n), m=m-1; end
    jj=1:m; n1j=(n+1-jj); bp(n1j)=-bp(jj); wf(n1j)=wf(jj);
    % end
    
 case 2, % p(x)=1/sqrt((x-a)*(b-x)), a=-1 and b=1 (default) 
  j=[1:n];
  wf = ones(1,n) * pi / n;
  bp=cos( (2*j-1)*pi / (2*n) );


 case 3, %p(x)=sqrt((x-a)*(b-x)),   a=-1   and b=1
  j=[1:n];
  wf = pi/ (n+1) *sin( j*pi / (n+1) ).^2;
  bp=cos( j*pi / (n+1) );


 case 4, %p(x)=sqrt((x-a)/(b-x)),   a=0   and b=1
    j=[1:n];
    bp=cos( (2*j-1)*pi /2/ (2*n+1) ).^2;
    wf=2*pi.*bp/(2*n+1) ;


 case 5, % %p(x)=1/sqrt(b-x),   a=0   and b=1
   [bp wf]=grule(2*n);
  wf(bp<0)=[];
  wf=wf*2;
   bp(bp<0)=[];
  bp=1-bp.^2;


 case 6, % %p(x)=sqrt(b-x),   a=0   and b=1
   [bp wf]=grule(2*n+1);
  wf(bp<=0)=[];
   bp(bp<=0)=[];
  wf=2*bp.^2.*wf;
  bp=1-bp.^2;
  
 case {7,8,9,10} ,%
  %7 p(x)=(x-a)^alpha*(b-x)^beta a=-1 b=1 Jacobi
  %8 p(x)=x^alpha*exp(-x) a=0,   b=inf generalized Laguerre
  %9 p(x)=exp(-x^2)       a=-inf, b=inf Hermite 
  %10 p(x)=1               a=-1 b=1        Legendre slower than 1
  % this procedure uses the coefficients a(j), b(j) of the
  %      recurrence relation
  %
  %           b p (x) = (x - a ) p   (x) - b   p   (x)
  %            j j            j   j-1       j-1 j-2
  %
  %      for the various classical (normalized) orthogonal polynomials,
  %      and the zero-th moment
  %
  %           muzero = integral w(x) dx
  %
  %      of the given polynomial's weight function w(x).  since the
  %      polynomials are orthonormalized, the tridiagonal matrix is
  %      guaranteed to be symmetric.
  %
  % 
  %         the input parameter alpha is used only for laguerre and
  %      jacobi polynomials, and the parameter beta is used only for
  %      jacobi polynomials.  the laguerre and jacobi polynomials
  %      require the gamma function.


  a=zeros(n,1);
  b=zeros(n-1,1);
  switch wfun
    case 7,  %jacobi
      ab = alpha + beta;
      abi = 2 + ab;
      muzero = 2^(ab + 1) * gamma(alpha + 1) * gamma(beta + 1) / gamma(abi);
      a(1) = (beta - alpha)/abi;
      b(1) = sqrt(4*(1 + alpha)*(1 + beta)/((abi + 1)*abi^2));
      a2b2 = beta^2 - alpha^2;
      
      i = (2:n-1)';
      abi = 2*i + ab;
      a(i) = a2b2./((abi - 2).*abi);
      a(n) =a2b2./((2*n - 2+ab).*(2*n+ab));
      b(i) = sqrt (4*i.*(i + alpha).*(i + beta)*(i + ab)./((abi.^2 - 1).*abi.^2));
   
    case 8, % Laguerre
      muzero=gamma(alpha+1);
      i = (1:n-1)';
      a(i) = 2 .* i - 1 + alpha;
      a(n)=2*n-1+alpha;
      b = sqrt( i .* (i + alpha) );
    case 9, %Hermite 
      i = (1:(n-1))';
      muzero = sqrt(pi);
      %a=zeros(m,1);
      b=sqrt(i/2);    
    case 10,  % legendre NB! much slower than wfun=1        
      muzero = 2;
      i = (1:n-1)';
      abi = i;
      b(i) = abi./sqrt(4*abi.^2 - 1);
      
  end
   
  %[v d] = eig( full(spdiags([b a b],-1:1,n,n )));
  [v d] = eig( diag(a) + diag(b,1) + diag(b,-1) );
  wf = v(1,:);
  if 1,
    [bp i] = sort( diag(d) );
    wf = wf(i);
  else % save some valuable time by not sorting
    bp = diag(d) ;
  end
  bp=bp';
  
  wf = muzero.* wf.^2;


otherwise, error('unknown weight function')
end

% end 