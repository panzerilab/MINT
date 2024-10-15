clc, clear all;

A = rand(2, 50,100);
B = rand(2, 50,100);

opts.tau = {[4,8],[5,9]};
opts.tpres = {35,41};

TE_v = TE({A,B}, opts);