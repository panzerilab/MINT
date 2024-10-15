clc, clear all, close all;
rng(0)
nTrials = 1000;
nTimepoints = 1000;
nNeurons = 4;

S1 = randi([0, 1], 1, nTrials, nTimepoints) * 2 - 1; 
S2 = randi([0, 1], 1, nTrials, nTimepoints) * 2 - 1; 


X = zeros(nNeurons, 5, nTrials, nTimepoints);

noiselevel = 0.5;
noiselevel_stim = 0.4;
Ex = noiselevel * randn(nNeurons, nTrials, nTimepoints);
Ex_stim = noiselevel_stim *  randn(nNeurons, nTrials, nTimepoints);

alpha = 1.2;
beta= 0.8;
c = 2.5;

X(:,5,:,:) = Ex;

X(1,1,:,1) = alpha*S1(:,1)+Ex_stim(1,:, 1);
X(2,2,:,1) =  beta*S2(:,1)+Ex_stim(2,:, 1);
X(4,1,:,1) = alpha*S1(:,1)+Ex_stim(4,:, 1);

for trial = 1:nTrials
    for tP = 2:nTimepoints
        X(1,1,trial,tP) = alpha*S1(1,trial,tP)+Ex_stim(1,trial, tP);
        X(2,2,trial,tP) = beta*S2(1,trial,tP)+Ex_stim(2,trial, tP);
        X(4,1,trial,tP) = alpha*S1(1,trial,tP)+Ex_stim(4,trial, tP);

        X(1,3,trial,tP) = c*X(2,2,trial,tP-1);
        X(2,3,trial,tP) = c*X(1,1,trial,tP-1);
        X(3,3,trial,tP) = c*X(1,1,trial,tP-1);
    end 
end 

X_resp = squeeze(sum(X, 2));

%% Information 

MI_opts.bin_methodX = 'eqpop';
MI_opts.bin_methodY = 'none';
MI_opts.n_binsX = 4;
MI_opts.method = 'dr';
MI_opts.bias = 'qe';
MI_opts.xtrp = 10;
MI_opts.verbose = false;

infoS1 = zeros(1,4);
infoS2 = zeros(1,4);

X_reshaped = reshape(X_resp, 4, nTrials*nTimepoints);
S1_reshaped = reshape(S1, 1, nTrials*nTimepoints);
S2_reshaped = reshape(S2, 1, nTrials*nTimepoints);

for i = 1:size(X_resp,1)
    infoS1(i) = cell2mat(information(X_reshaped(i,:), S1_reshaped, MI_opts,{'I'}));
    infoS2(i) = cell2mat(information(X_reshaped(i,:), S2_reshaped, MI_opts,{'I'}));
end 

figure;
subplot(1,2,1);
bar(infoS1);
title('Mutual Information between Channels and S1');
xlabel('Channels');
ylabel('Mutual Information');
set(gca, 'XTick', 1:4, 'XTickLabel', {'X1', 'X2', 'X3', 'X4'});

subplot(1,2,2);
bar(infoS2);
title('Mutual Information between Channels and S2');
xlabel('Channels');
ylabel('Mutual Information');
set(gca, 'XTick', 1:4, 'XTickLabel', {'X1', 'X2', 'X3', 'X4'});

%% TE
TE = zeros(4,4);

TE_opts.bin_methodX = 'eqpop';
TE_opts.bin_methodY = 'eqpop';
TE_opts.n_binsX = 4;
TE_opts.n_binsY = 4;
TE_opts.method = 'dr';
TE_opts.bias = 'pt';
TE_opts.xtrp = 10;
TE_opts.verbose = false;
TE_opts.taux = -1;
TE_opts.tauy = -1;

X_perm = permute(X_resp, [1, 3, 2]);


for x1 = 1:size(X_resp,1)
    for x2 = 1:size(X_resp,1)
        if x1 ~= x2
             TE(x1, x2) = cell2mat(transferentropy(squeeze(X_perm(x1,:,:)),squeeze(X_perm(x2,:,:)),TE_opts, {'TE'})); 
        end 
    end 
end 

%% FIT

FIT_opts.bin_methodX = 'eqpop';
FIT_opts.bin_methodY = 'eqpop';
FIT_opts.bin_methodS = 'none';
FIT_opts.n_binsX = 4;
FIT_opts.n_binsY = 4;
FIT_opts.method = 'dr';
FIT_opts.bias = 'naive';
%FIT_opts.xtrp = 10;
FIT_opts.verbose = false;

FIT1 = zeros(4,4);
FIT2 = zeros(4,4);

for x1 = 1:size(X_resp,1)
    for x2 = 1:size(X_resp,1)
        hX = reshape(X_resp(x1, :, 1:end-1), 1, (nTimepoints-1)*nTrials);
        hY = reshape(X_resp(x2, :, 1:end-1), 1, (nTimepoints-1)*nTrials);
        Y  = reshape(X_resp(x2, :, 2:end), 1, (nTimepoints-1)*nTrials);
        S1_FIT = reshape(S1(1,:,2:end), 1, (nTimepoints-1)*nTrials);
        S2_FIT = reshape(S2(1,:,2:end), 1, (nTimepoints-1)*nTrials);
        if x1 ~= x2
             FIT1(x1, x2) = FIT(S1_FIT, hX, hY, Y, FIT_opts); 
             FIT2(x1, x2) = FIT(S2_FIT, hX, hY, Y, FIT_opts); 
        end 
    end 
end 

