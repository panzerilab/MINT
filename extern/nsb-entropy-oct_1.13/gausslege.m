function I=gausslege(f,a,b,n)

%I=gaussleg(f,a,b,n)
%Aproximates integral using Gauss-Legendre cuadrature method

%Legendre polin
p=polegende(n);
%Polin roots
x=roots(p(n+1,:));

%Change of integration variable if it's needed
if a~=-1 | b~=1
   y=flege(f,a,b);
   G=subs(y,x);
else
   G=feval(f,x);		
end


%Polin derivate
pn=polyder(p(n+1,:));

%Calculus of the coeficients
for i=1:n
   C(i)=2./((1-x(i).^2).*((polyval(pn,x(i))).^2));
end

%Escalar product of the function at the nodes and the coeficients
I=dot(C,G);


