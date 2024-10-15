function [II_v, II_v_uncorrected] = iispeed()
rng('default');
S = rand(1, 100);
R = rand(1, 100);
C = rand(1, 100);
opts.bias = 'qe';
opts.parallel = 1;
opts.xtrp= 50;
opts.bin_methodS = 'eqpop';
opts.bin_methodR = 'eqpop';
opts.bin_methodC = 'eqpop';
opts.n_binsS = 2;
opts.n_binsR = 2;
opts.n_binsC = 2;

[II_v, II_v_uncorrected] = II(S,R,C, opts);
% [I_v, PID_SR_C, PID_RC_S] =  II(S,R,C, opts);
end