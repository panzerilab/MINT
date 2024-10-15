clc, clear,

% -------------------------------------------------------------------------
%                            Tutorial MINT
% -------------------------------------------------------------------------
% In this tutorial, 

% -------------------------------------------------------------------------
% Step 1: Initialize Parameters and Options
% -------------------------------------------------------------------------

rng(0);
 
nNeurons = 100;
xi = linspace(0, 180, nNeurons);
sigma = 10;
s_d = 0:1:180;
nTimepoints = 50;

response_d = zeros(nNeurons, length(s_d));
for i = 1:length(s_d)
     response_d(:, i) = exp(-((s_d(i) - xi).^2) / (2 * sigma^2));
end 

s_i = 0:0.05:1;
response_i = zeros(nNeurons, length(s_i));
iota = linspace(0, 1, nNeurons);
sigma = 0.1;
for i = 1:length(s_i)
     response_i(:, i) = exp(-((s_i(i) - iota).^2) / (2 * sigma^2));
end

w0 = 2;  
p = 0.1;
sigma_w = 50;

connections = rand(nNeurons, nNeurons) < p;
connections(eye(size(connections)) == 1) = 0;
wij = zeros(nNeurons, nNeurons);

for i = 1:nNeurons
    for j = 1:nNeurons
        deltaXi = xi(i) - xi(j);
        wij(i, j) = w0 * exp(-(deltaXi.^2) / (2 * sigma_w^2));
    end
end
wij = wij .* connections;
v_d = wij*response_d;
v_i = wij*response_i;


s = [randi([0, 180],1, 150); rand(1,150)];

for neuron = 1:nNeurons
     for t = 1:length(s)
         intensity= s(2,t);
         diff = abs(s_i - intensity);
        [~, index] = min(diff);
         direction = s(1,t);
        firing_rate = v_d(neuron, direction) + v_i(neuron, index);
        spikes = poisson_spike_gen(1:nTimepoints, firing_rate, 0.1);
        R(neuron,:, t) = spikes;
     end 
end 

for t = 1:length(s) 
    for nN = 1:nNeurons
        R_spike_sum(nN, t) = sum(R(nN, :, t));
    end
end

R = squeeze(sum(R, 2)/nTimepoints);
opts.bias = 'naive';
opts.method = 'dr';
outputsList = {'HX', 'HXY'};
opts.bin_methodX = 'eqspace';
opts.n_binsX = 5;
opts.bin_methodY = 'eqspace';
opts.n_binsY = 2;
opts.verbose = 0;
opts.tauy = -0.5;
opts.taux = -0.5;

I_d = zeros(1, nNeurons);
for n = 1:nNeurons
    outputs = information(R(n,:), s(1,:), opts, {'I'});
    I_d(n) = outputs{1};
end 

I_i = zeros(1, nNeurons);
for n = 1:nNeurons
    outputs = information(R(n,:), s(2,:), opts, {'I'});
    I_i(n) = outputs{1};
end 
info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
pairIdx = 0;
for cell1 = 1
        for cell2 = 4
            pairIdx = pairIdx+1;
            jointResp = [R_spike_sum(cell1,:);R_spike_sum(cell2,:)];
            infoBdw = information(jointResp,s,opts,{'I','ILIN','ISS','ICI','ICD'});
            for bdwIdx = 1:numel(info_bdw_terms)
                bdwLab = info_bdw_terms{bdwIdx};
                MI_breakdown.(bdwLab)(pairIdx) = infoBdw{bdwIdx};
            end
        end
end

for cell1 = 1
        for cell2 = 3
            pairIdx = pairIdx+1;
            jointResp = [R_spike_sum(cell1,:);R_spike_sum(cell2,:)];
            infoBdw = information(jointResp,s,opts,{'I','ILIN','ISS','ICI','ICD'});
            for bdwIdx = 1:numel(info_bdw_terms)
                bdwLab = info_bdw_terms{bdwIdx};
                MI_breakdown.(bdwLab)(pairIdx) = infoBdw{bdwIdx};
            end
        end
end
%%
[meanJoint, errJoint] = compute_stats(MI_breakdown.Joint);
[meanLin, errLin] = compute_stats(MI_breakdown.ILIN);
[meanISS, errISS] = compute_stats(MI_breakdown.ISS);
[meanICI, errICI] = compute_stats(MI_breakdown.ICI);
[meanICD, errICD] = compute_stats(MI_breakdown.ICD);


fSize = 15;
figure('Position', [360, 347, 785, 271])
hold on 
% Barplots
% Joint information
h{1} = bar(1, meanJoint,'FaceColor',[0.4 0.4 0.4]);
errorbar(1, meanJoint,errJoint,'k--')
% I_Lin
h{2} = bar(2, meanLin,'FaceColor',[0.8 0.4 0.4]);
errorbar(2, meanLin,errLin,'k--')
% Sig. sim.
h{3} = bar(3, meanISS,'FaceColor',[0.4 0.8 0.4]);
errorbar(3, meanISS,errISS,'k--')
% Cor-indep
h{4} = bar(4, meanICI,'FaceColor',[0.4 0.4 0.8]);
errorbar(4, meanICI,errICI,'k--');
% Cor-dep
h{5} = bar(5, meanICD,'FaceColor',[0.8 0.4 0.8]);
errorbar(5, meanICD,errICD,'k--');

xticks([]);
set(gca, 'FontSize', 12);
ylabel('Information [bits]')
legend_obj = legend([h{1} h{2} h{3} h{4} h{5}], 'Joint', 'Linear', 'Sig. sim.', 'Corr. indep.', 'Corr. dep.');
legend_obj.FontSize = 14;
legend('AutoUpdate','off');
ylim([-0.1, 0.7])
set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 60]);
ax = gca;
ax.FontSize = 12;
title('Information Breakdown', 'FontSize', fSize+5);



function [mean_val, sem_val] = compute_stats(data)
    mean_val = mean(data);
    std_val = std(data);
    n = length(data);
    sem_val = std_val / sqrt(n);
end

