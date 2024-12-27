n = [4 2 3 0 2 4 0 0 2]; 
K = length(n);
nx = n(n>0);
kx = ones(size(nx));
qfun=1
precision=.1
todo=1
[S_nsb, dS_nsb, S_cl, dS_cl, xi_cl, S_ml,errcode]=find_nsb_entropy (kx, nx, K, precision,qfun)
