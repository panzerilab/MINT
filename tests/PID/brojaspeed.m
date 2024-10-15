function [pid_v] = brojaspeed()
% rng('default');
dimsize = 8;
prob_xyz = rand([dimsize, dimsize, dimsize]);
prob_xyz = prob_xyz /sum(prob_xyz,'all');

pid_v = pidBROJA(prob_xyz);
% [PID_v] = PID(Y, X1, X2, opts);
end