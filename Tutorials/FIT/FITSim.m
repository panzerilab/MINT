clear all, close all;
rng(0);


simReps = 3;

nNeurons = 4;
nTrials = 26;
nTimepoints = 1000;
delay = 6;

alpha = 1;
c = 1;

X = zeros(nNeurons, 2, nTimepoints, nTrials);
infoS1 = zeros(4,nTimepoints, simReps);
infoS2 = zeros(4,nTimepoints, simReps);

FIT1 = zeros(4,4, nTrials, simReps);
FIT2 = zeros(4,4, nTrials, simReps);

TE = zeros(4,4, simReps);

opts = [];
opts.verbose = false;
opts.method = "dr";
opts.bias = 'naive';
opts.xtrp = 3;
opts.btsp = 0;
opts.n_binsS = 2; % Number of stimulus values
opts.n_binsX = 2;
opts.n_binsY = 2;
opts.bin_methodS = 'none';
opts.bin_methodX = 'eqpop';
opts.bin_methodY = 'eqpop';
opts.taux = [-delay, -delay+1];
opts.tauy = -1;
opts.btsp_variables = {};
opts.parallel = 0;

for repIdx = 1:simReps
    repIdx
    X = zeros(nNeurons, 2, nTimepoints, nTrials);
    S1 = 2*randi([0, 1], nTimepoints, nTrials) - 1;
    S2 = 2*randi([0, 1], nTimepoints, nTrials) - 1;

    epsX = 0.5;   

    for t=1:nTimepoints

        X(:,:,t,:) = epsX*randn(4,2,1, nTrials);
        
        X(1,1,t,:) = squeeze(X(1,1,t,:)) + (alpha*S1(t,:))';
        X(2,2,t,:) = squeeze(X(2,2,t,:)) + (alpha*S2(t,:))';
        X(4,1,t,:) = squeeze(X(4,1,t,:)) + (alpha*S1(t,:))';
       
 
        if t > delay 
            X(1,2,t,:) = X(1,2,t,:) + c*X(2,2,t-delay,:);
            X(2,1,t,:) = X(2,1,t,:) + c*X(1,1,t-delay,:);
            X(3,1,t,:) = X(3,1,t,:) + c*X(1,1,t-delay,:);
            X(4,2,t,:) = X(4,2,t,:) + 3*c*X(3,2,t-delay,:);
        end
        
        opts.bin_methodY = 'none';

        for neuron = 1:nNeurons
            infoS1(neuron, t, repIdx) = cell2mat(information(squeeze(sum(X(neuron,:,t,:),2)), S1(t,:) ,opts,{'I'}));
            infoS2(neuron, t, repIdx) = cell2mat(information(squeeze(sumX(neuron,:,t,:)), S2(t,:) ,opts,{'I'}));
        end

    end

    opts.bin_methodY = 'eqpop';
    opts.bias = 'naive';
    for neuron1 = 1:nNeurons
        for neuron2 = 1:nNeurons
            if neuron1 ~= neuron2
                TE(neuron1, neuron2, repIdx) = cell2mat(transferentropy(squeeze(sum(X(neuron1,:,:,:),2)), squeeze(sum(X(neuron2,:,:,:),2)),opts, {'TE'}));          
            end
        end
    end

    for trial = 1:nTrials   

        X_pres = X(:,:,delay+1:end,trial); 
        X_past = X(:,:,1:end-delay,trial); 
    
        opts.bias = 'naive';
         for neuron1 = 1:nNeurons
            for neuron2 = 1:nNeurons
                if neuron1 ~= neuron2
                     FIT1(neuron1, neuron2, trial, repIdx) = FIT(S1(delay+1:end,trial)',squeeze(sum(X_past(neuron1,:,:,:),2))',squeeze(sum(X_past(neuron2,:,:,:),2))',squeeze(sum(X_pres(neuron2,:,:,:),2))',opts);            
                     FIT2(neuron1, neuron2, trial, repIdx) = FIT(S2(delay+1:end,trial)',squeeze(sum(X_past(neuron1,:,:,:),2))',squeeze(sum(X_past(neuron2,:,:,:),2))',squeeze(sum(X_pres(neuron2,:,:,:),2))',opts);            
                end 
            end 
        end 
    end 
end
%%
infoS1_mean = mean(mean(infoS1,2),3);
infoS2_mean = mean(mean(infoS2,2),3);
FIT1_mean = mean(mean(FIT1, 3),4);
FIT2_mean = mean(mean(FIT2, 3),4);

TE_mean = mean(TE, 3);

%% Plot 
X_labels = {'X1', 'X2', 'X3', 'X4'};
data = [infoS1_mean(:) infoS2_mean(:)];

figure;
b = bar(data);
b(1).FaceColor = [0.4 0.4 0.4];
b(2).FaceColor = [0.8 0.4 0.4];

set(gca, 'XTickLabel', X_labels);
xlabel('');
ylabel('Mean Information');
legend({'I(S1;X)','I(S2;X)'}, 'Location', 'northeast');
title('Info about S1 and S2');
saveas(gcf, 'mean_information_plot.png');

figure 
heatmap(TE_mean);
xlabel('Target Neuron');
ylabel('Source Neuron');
title('Transfer Entropy (TE)');
saveas(gcf, 'Heatmap_TE.png');

figure 
heatmap(FIT1_mean);
xlabel('Target Neuron');
ylabel('Source Neuron');
title('FIT S1');
saveas(gcf, 'HEatmap_FIT1.png');

figure 
heatmap(FIT2_mean);
xlabel('Target Neuron');
ylabel('Source Neuron');
title('FIT S2');
saveas(gcf, 'HEatmap_FIT2.png');
