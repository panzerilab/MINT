function [f] = psi(z)
%%Psi     Psi (or Digamma) function valid in the entire complex plane.
%%
%%                 d
%%        Psi(z) = --log(Gamma(z))
%%                 dz
%%
%%usage: [f] = psi(z)
%%
%%tested under version 5.3.1
%%
%%        Z may be complex and of any size.
%%
%%        This program uses the analytical derivative of the
%%        Log of an excellent Lanczos series approximation
%%        for the Gamma function.
%%        
%%References: C. Lanczos, SIAM JNA  1, 1964. pp. 86-96
%%            Y. Luke, "The Special ... approximations", 1969 pp. 29-31
%%            Y. Luke, "Algorithms ... functions", 1977
%%            J. Spouge,  SIAM JNA 31, 1994. pp. 931
%%            W. Press,  "Numerical Recipes"
%%            S. Chang, "Computation of special functions", 1996
%%
%%
%%see also:   GAMMA GAMMALN GAMMAINC
%%see also:   mhelp psi
%%see also:   mhelp GAMMA


%%Paul Godfrey
%%pgodfrey@intersil.com
%%11-15-00
% 

 [row, col] = size(z);
 z=z(:);
 zz=z;


f = 0.*z; %# reserve space in advance


%reflection point
p=find(real(z)<0.5);
if ~isempty(p)
   z(p)=1-z(p);
end;


%##Lanczos approximation for the complex plane
c = [   .9999999999998099322768;
     676.5203681218850985670;
   -1259.139216722402870472;
     771.3234287776530788487;
    -176.6150291621405990658;
      12.50734327868690481446;
      -0.1385710952657201168955;
       0.9984369578019570859563e-5;
       0.1505632735149311558338e-6];
 
g=7;
%##see gamma for calculation details...


n=0;
d=0;
for k=g+2:-1:2
    dz=1./(z+k-2);
    dd=c(k).*dz;
    d=d+dd;
    n=n-dd.*dz;
end;
d=d+c(1);
gg=z+g-0.5;
f = n./d + (z-0.5)./gg - 1 + log(gg);


if ~isempty(p)
   f(p) = f(p)-pi*cot(pi*zz(p));
end;


p=find(round(zz)==zz & real(zz)<=0 & imag(zz)==0);
if ~isempty(p)
   f(p) = Inf;
end;


f=reshape(f,row,col);

