clear all;
rng(0);

tic;
simReps = 5; % number of repetitions of the simulation

nNeurons = 4;
nTrials = 200;
nTimepoints = 30;
stim_window = [3 8];
delay = 7;
nShuff = 3; % number of shufflings
null_samples = 200; % number of null hypothesis samples obtained by bootstrapping shufflings across repetitions

c = 1;

X = zeros(nNeurons, 2, nTimepoints, nTrials);
infoS1 = zeros(4,nTimepoints, simReps);
infoS2 = zeros(4,nTimepoints, simReps);

FIT1 = zeros(4,4, nTimepoints, simReps);
FIT2 = zeros(4,4, nTimepoints, simReps);
FIT1_sh.conditioned = zeros(4,4, nTimepoints, simReps, nShuff);
FIT2_sh.conditioned = zeros(4,4, nTimepoints, simReps, nShuff);
FIT1_sh.simple = zeros(4,4, nTimepoints, simReps, nShuff);
FIT2_sh.simple = zeros(4,4, nTimepoints, simReps, nShuff);

TE = zeros(4,4, nTimepoints, simReps);

opts = [];
opts.verbose = false;
opts.method = "dr";
opts.bias = 'qe';
opts.xtrp = 3;
opts.btsp = 0;
opts.n_binsS = 2; % Number of stimulus values
opts.n_binsX = 2;
opts.n_binsY = 2;
opts.bin_methodS = 'none';
opts.bin_methodX = 'eqpop';
opts.bin_methodY = 'eqpop';
opts.taux = -delay;
opts.tauy = -1;
opts.btsp_variables = {};
opts.parallel = 0;

for repIdx = 1:simReps
    repIdx
    X = zeros(nNeurons, 2, nTimepoints, nTrials);
    S1 = randi(opts.n_binsS,1,nTrials);
    S2 = randi(opts.n_binsS,1,nTrials);

    epsX = 0.5;

    for t=1:nTimepoints

        X(:,:,t,:) = epsX*randn(4,2,1, nTrials);


        % if (t > stim_window(1) || t > stim_window(3))  && (t <  stim_window(2) || t < stim_window(4))
        if (t > stim_window(1))  && (t <  stim_window(2))
            X(1,1,t,:) = squeeze(X(1,1,t,:)) +epsX*randn+2*(S1-1.5)';
            X(2,2,t,:) = squeeze(X(2,2,t,:)) +epsX*randn+2*(S2-1.5)';
            X(4,1,t,:) = squeeze(X(4,1,t,:)) +epsX*randn+2*(S1-1.5)';
        end

        if t > delay
            X(1,2,t,:) = X(1,2,t,:) + c*X(2,2,t-delay,:);
            X(2,1,t,:) = X(2,1,t,:) + c*X(1,1,t-delay,:);
            X(3,1,t,:) = X(3,1,t,:) + c*X(1,1,t-delay,:);
            X(4,2,t,:) = X(4,2,t,:) + 3*c*X(3,2,t-delay,:); % noise transfer
        end

        opts.bin_methodY = 'none';

        for neuron = 1:nNeurons
            infoS1(neuron, t, repIdx) = cell2mat(information(squeeze(sum(X(neuron,:,t,:),2))', S1 ,opts,{'I'}));
            infoS2(neuron, t, repIdx) = cell2mat(information(squeeze(sum(X(neuron,:,t,:),2))', S2 ,opts,{'I'}));
        end

    end

    for t = delay+1:nTimepoints
        opts.bin_methodY = 'eqpop';

        X_pres = X(:,:,t,:);
        X_past = X(:,:,t-delay,:);

        for neuron1 = 1:nNeurons
            for neuron2 = 1:nNeurons
                if neuron1 ~= neuron2
                    opts.bias = 'qe';
                    FIT1(neuron1, neuron2,t, repIdx) = FIT(S1,squeeze(sum(X_past(neuron1,:,:,:),2))',squeeze(sum(X_past(neuron2,:,:,:),2))',squeeze(sum(X_pres(neuron2,:,:,:),2))',opts);
                    FIT2(neuron1, neuron2,t, repIdx) = FIT(S2,squeeze(sum(X_past(neuron1,:,:,:),2))',squeeze(sum(X_past(neuron2,:,:,:),2))',squeeze(sum(X_pres(neuron2,:,:,:),2))',opts);
                    opts.bias = 'naive';
                    TE(neuron1, neuron2,t, repIdx) = cell2mat(transferentropy(squeeze(sum(X(neuron1,:,:,:),2)), squeeze(sum(X(neuron2,:,:,:),2)),opts, {'TE'}));
                end
            end
        end


        % Compute null hypothesis for FIT and TE
        for shIdx = 1:nShuff

            shuffledS1 = S1(randperm(numel(S1)));
            shuffledS2 = S2(randperm(numel(S2)));

            for neuron1 = 1:nNeurons

                % Permute data to have trials in the first dimension (as
                % required by the shuffle function)
                sender = X(neuron1,:,:,:);
                sender_perm = permute(sender,[4 2 3 1]);

                % Unconstrained shuffling of the sender
                sender_all_shuff_perm = shuffle(S1',sender_perm,0,[1 0 0 0]);
                sender_all_shuff = permute(sender_all_shuff_perm,[4 2 3 1]);

                % Shuffle the sender at fixed values of S1
                sender_s1_shuff_perm = shuffle(S1',sender_perm,1,[1 0 0 0]);
                sender_s1_shuff = permute(sender_s1_shuff_perm,[4 2 3 1]);

                % Shuffle the sender at fixed values of S2
                sender_s2_shuff_perm = shuffle(S2',sender_perm,1,[1 0 0 0]);
                sender_s2_shuff = permute(sender_s2_shuff_perm,[4 2 3 1]);

                X_s1_shuff = X; X_s2_shuff = X; X_all_shuff = X;
                X_s1_shuff(neuron1,:,:,:) = sender_s1_shuff;
                X_s2_shuff(neuron1,:,:,:) = sender_s2_shuff;
                X_all_shuff(neuron1,:,:,:) = sender_all_shuff;

                % compute past and present of neural signals
                X_all_shuff_pres = X_all_shuff(:,:,t,:);
                X_all_shuff_past = X_all_shuff(:,:,t-delay,:);
                X_s1_shuff_pres = X_s1_shuff(:,:,t,:);
                X_s1_shuff_past = X_s1_shuff(:,:,t-delay,:);
                X_s2_shuff_pres = X_s2_shuff(:,:,t,:);
                X_s2_shuff_past = X_s2_shuff(:,:,t-delay,:);

                for neuron2 = 1:nNeurons
                    if neuron1 ~= neuron2
                        opts.bias = 'qe';
                        % Conditioned shuffling (fixed stimulus)
                        FIT1_sh.conditioned(neuron1, neuron2,t, repIdx, shIdx) = FIT(S1,squeeze(sum(X_s1_shuff_past(neuron1,:,:,:),2))',squeeze(sum(X_s1_shuff_past(neuron2,:,:,:),2))',squeeze(sum(X_s1_shuff_pres(neuron2,:,:,:),2))',opts);
                        FIT2_sh.conditioned(neuron1, neuron2,t, repIdx, shIdx) = FIT(S2,squeeze(sum(X_s2_shuff_past(neuron1,:,:,:),2))',squeeze(sum(X_s2_shuff_past(neuron2,:,:,:),2))',squeeze(sum(X_s2_shuff_pres(neuron2,:,:,:),2))',opts);

                        % Simple S shuffling
                        FIT1_sh.simple(neuron1, neuron2,t, repIdx, shIdx) = FIT(shuffledS1,squeeze(sum(X_past(neuron1,:,:,:),2))',squeeze(sum(X_past(neuron2,:,:,:),2))',squeeze(sum(X_pres(neuron2,:,:,:),2))',opts);
                        FIT2_sh.simple(neuron1, neuron2,t, repIdx, shIdx) = FIT(shuffledS2,squeeze(sum(X_past(neuron1,:,:,:),2))',squeeze(sum(X_past(neuron2,:,:,:),2))',squeeze(sum(X_pres(neuron2,:,:,:),2))',opts);
                        opts.bias = 'naive';
                        TE_sh(neuron1, neuron2, repIdx,t, shIdx) = cell2mat(transferentropy(squeeze(sum(X_all_shuff(neuron1,:,:,:),2)), squeeze(sum(X_all_shuff(neuron2,:,:,:),2)),opts, {'TE'}));                
                    end
                end
            end
        end
    end
end
%%
infoS1_mean = mean(mean(infoS1,2),3);
infoS2_mean = mean(mean(infoS2,2),3);
FIT1_mean = mean(mean(FIT1,3),4);
FIT2_mean = mean(mean(FIT2,3),4);
TE_mean = mean(mean(TE, 3),4);
%%

% Bootstrap shuffled values across repetitions of the simulations to boost
% the number of null hypothesis samples
FIT1_sh_perm.simple = permute(FIT1_sh.simple,[4 1 2 3 5]);
FIT1_sh_mean.simple = btstrp_shuff(FIT1_sh_perm.simple,null_samples);
FIT2_sh_perm.simple = permute(FIT2_sh.simple,[4 1 2 3 5]);
FIT2_sh_mean.simple = btstrp_shuff(FIT2_sh_perm.simple,null_samples);

FIT1_sh_perm.conditioned = permute(FIT1_sh.conditioned,[4 1 2 3 5]);
FIT1_sh_mean.conditioned = btstrp_shuff(FIT1_sh_perm.conditioned,null_samples);
FIT2_sh_perm.conditioned = permute(FIT2_sh.conditioned,[4 1 2 3 5]);
FIT2_sh_mean.conditioned = btstrp_shuff(FIT2_sh_perm.conditioned,null_samples);

% Take the element-wise maximum across the two null hypotheses
FIT1_sh_mean.max = max(FIT1_sh_mean.simple,FIT1_sh_mean.conditioned);
FIT2_sh_mean.max = max(FIT2_sh_mean.simple,FIT2_sh_mean.conditioned);

TE_sh_perm = permute(TE_sh,[4 1 2 3 5]);
TE_sh_mean = btstrp_shuff(TE_sh_perm,100);

% Determine significance thresholds for FIT 1 and FIT 2
sig_thresh_FIT1 = prctile(FIT1_sh_mean.max,95,3);
sig_thresh_FIT2 = prctile(FIT2_sh_mean.max,95,3);
sig_thresh_TE = prctile(TE_sh_mean,95,3);

toc

save('FIT_network_results_120624')
%%
% Info encoded in each neurons
figure('Position',[270,197,798,277])
subplot(1,2,1)
plot(squeeze(mean(infoS1,3))','LineWidth',2)
legend('X1','X2','X3','X4','AutoUpdate','off')
xlabel('time')
ylabel('[bits]')
title('info S1')

subplot(1,2,2)
plot(squeeze(mean(infoS2,3))','LineWidth',2)
xlabel('time')
ylabel('[bits]')
title('info S2')

saveas(gcf, 'Information_FITsim.png')

% Plot TE and FIT
figure('Position',[25,230,1179,296])
subplot(1,3,1)
imagesc(TE_mean)
for i = 1:nNeurons
    for j = 1:nNeurons
        if(TE_mean(i,j)>sig_thresh_TE(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
yticks(1:4)
yticklabels({'1','2','3','4'})
title('TE')
colorbar()

subplot(1,3,2)
imagesc(FIT1_mean)
for i = 1:nNeurons
    for j = 1:nNeurons
        if(FIT1_mean(i,j)>sig_thresh_FIT1(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
title('FIT_{S1}')
yticks(1:4)
yticklabels({'1','2','3','4'})
colorbar()

subplot(1,3,3)
imagesc(FIT2_mean)
for i = 1:nNeurons
    for j = 1:nNeurons
        if(FIT2_mean(i,j)>sig_thresh_FIT2(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
yticks(1:4)
yticklabels({'1','2','3','4'})
title('FIT_{S2}')
colorbar()

saveas(gcf, 'TE_and_FIT_FITsim.png')