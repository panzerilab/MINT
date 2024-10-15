% Bidirectional transmission of info about different features

clear all; %close all;

rng(0) % For reproducibility
save_results = 0; % set to 1 to save results file
% results_path =  'C:\Users\mcelotto\Desktop\Neural Computation\PhD\Various\Roberto\NIT_paper\Fig4_FIT'; % path to the directory to save results
% MINT_path = 'C:\Users\mcelotto\Desktop\Neural Computation\PhD\Sur_project\scripts\github_sync\MINT';
% addpath(MINT_path)

% Simulation parameters
nTrials_per_stim = [1000]; % number of trials per stimulus value (power)
simReps = 50; % repetitions of the simulation
epsY = 1; 
epsX = 1; % standard deviation of bagaussian noise in X and Y

tparams.simLen = 40; % simulation time, in units of 10ms
tparams.stimWin = [10 15]; % X stimulus encoding window, in units of 10ms
tparams.delay = [7]; % communication delays, in units of 10ms

% Define information options
opts = [];
opts.verbose = false;
opts.method = "dr";
opts.bias = 'naive';
opts.btsp = 0;
opts.n_binsS = 2; % Number of stimulus values
opts.n_binsX = 2;
opts.n_binsY = 2;
opts.bin_methodS = 'none';
opts.bin_methodX = 'eqpop';
opts.bin_methodY = 'eqpop';
opts.btsp_variables = {};

% Connectivity weights
w_xy = 1; % strength of X to Y transfer
w_yx = 1; % strength of Y to X transfer

% Initialize structures
fit.X2Y.S1 = nan(1,simReps); fit.X2Y.S2 = fit.X2Y.S1; te.X2Y = fit.X2Y.S1; 
fit.Y2X.S1 = nan(1,simReps); fit.Y2X.S2 = fit.Y2X.S1; te.Y2X = fit.Y2X.S1; 
infohX.S1 = te.X2Y; infohY.S1 = te.X2Y; infoYt.S1 = te.X2Y; infoXt.S1 = te.X2Y;
infohX.S2 = te.X2Y; infohY.S2 = te.X2Y; infoYt.S2 = te.X2Y; infoXt.S2 = te.X2Y;
%% Run simulation

nTrials = (nTrials_per_stim)*opts.n_binsS; % Compute number of trials

for repIdx = 1:simReps
    repIdx

    % Draw the stimulus value for each trial
    S1 = randi(opts.n_binsS,1,nTrials);
    S2 = randi(opts.n_binsS,1,nTrials);

    for t=1:tparams.simLen

        % sBaseline noise
        X_S1(t,:) = epsX*randn(1,nTrials); % X noise time series
        Y_S1(t,:) = epsY*randn(1,nTrials); % X noise time series
        X_S2(t,:) = epsX*randn(1,nTrials); % X noise time series
        Y_S2(t,:) = epsY*randn(1,nTrials); % X noise time series

        % Stimulus encoding
        if t > tparams.stimWin(1) & t<tparams.stimWin(2)
            X_S1(t,:) = X_S1(t,:) + 2*(S1-1.5);
            Y_S2(t,:) = Y_S2(t,:) + 2*(S2-1.5);
        end
    
        % Communication between the two areas
        if t > tparams.delay 
            % Transfer of S1 from X to Y (in the S1 dimension)
            Y_S1(t,:) = Y_S1(t,:) + w_xy*X_S1(t-tparams.delay,:);
            % Transfer of S2 from Y to X (in the S2 dimension)
            X_S2(t,:) = X_S2(t,:) + w_yx*Y_S2(t-tparams.delay,:);
        end

        opts.bin_methodX = 'none'; % here X is the stimulus (already binary)

        infoX.S1(repIdx,t) = cell2mat(information(S1,[X_S1(t,:);X_S2(t,:)],opts,{'I'}));
        infoX.S2(repIdx,t) = cell2mat(information(S2,[X_S1(t,:);X_S2(t,:)],opts,{'I'}));
        infoY.S1(repIdx,t) = cell2mat(information(S1,[Y_S1(t,:);Y_S2(t,:)],opts,{'I'}));
        infoY.S2(repIdx,t) = cell2mat(information(S2,[Y_S1(t,:);Y_S2(t,:)],opts,{'I'}));

    end

    % First time point at which Y receives stim info from X
    t = tparams.stimWin(1)+tparams.delay+1; % first encoding time point (t = 150ms) + delay (70ms)
    d = tparams.delay;

    opts.bin_methodX = 'eqpop';

    % Define 2D signals for the past/present of X/Y combining the two dimensions
    joint_X = cat(3,X_S1,X_S2);
    joint_Y = cat(3,Y_S1,Y_S2);
    joint_X=permute(joint_X,[3,1,2]);
    joint_Y=permute(joint_Y,[3,1,2]);

    joint_X_pres = squeeze(joint_X(:,t,:));
    joint_Y_pres = squeeze(joint_Y(:,t,:));
    joint_X_past = squeeze(joint_X(:,t-d,:));
    joint_Y_past = squeeze(joint_Y(:,t-d,:));

    % Compute FIT and TE
    fit.X2Y.S1(repIdx)= FIT(S1,joint_X_past,joint_Y_past,joint_Y_pres,opts);
    fit.X2Y.S2(repIdx)= FIT(S2,joint_X_past,joint_Y_past,joint_Y_pres,opts);
    fit.Y2X.S1(repIdx)= FIT(S1,joint_Y_past,joint_X_past,joint_X_pres,opts);
    fit.Y2X.S2(repIdx)= FIT(S2,joint_Y_past,joint_X_past,joint_X_pres,opts);

    te.X2Y(repIdx)= cell2mat(information([joint_X_past;joint_Y_past],joint_Y_pres,opts,{'I'}))-cell2mat(information(joint_Y_past,joint_Y_pres,opts,{'I'}));
    te.Y2X(repIdx)= cell2mat(information([joint_X_past;joint_Y_past],joint_X_pres,opts,{'I'}))-cell2mat(information(joint_X_past,joint_X_pres,opts,{'I'}));
end

if save_results
    fname = ['\FIT_simulation_50reps_100624.mat'];
    save([results_path,fname])
end

%% Plots

% Plot info profiles
figure()
subplot(2,1,1)
hold on
plot(1:tparams.simLen,mean(infoX.S1),'linewidth',2)
plot(1:tparams.simLen,mean(infoY.S1),'linewidth',2)
ylabel('[bits]')
xlabel('time')
title('Info about S1')
legend('I(S;X)','I(S;Y)')

subplot(2,1,2)
hold on
plot(1:tparams.simLen,mean(infoX.S2),'linewidth',2)
plot(1:tparams.simLen,mean(infoY.S2),'linewidth',2)
ylabel('[bits]')
xlabel('time')
title('Info about S2')

% TE figure
figure()
hold on
bar(1,mean(te.X2Y))
errorbar(1,mean(te.X2Y),std(te.X2Y)/sqrt(numel(simReps)),'ko')
bar(2,mean(te.Y2X))
errorbar(2,mean(te.Y2X),std(te.Y2X)/sqrt(numel(simReps)),'ko')

set(gca, 'XTick',[1,2], 'XTickLabel', {'X->Y','Y->X'})
ylabel('[bits]')
title('TE in the two directions')

% FIT figure
figure()
hold on
bar(1,mean(fit.X2Y.S1))
errorbar(1,mean(fit.X2Y.S1),std(fit.X2Y.S1)/sqrt(numel(simReps)),'ko')
bar(2,mean(fit.Y2X.S1))
errorbar(2,mean(fit.Y2X.S1),std(fit.Y2X.S1)/sqrt(numel(simReps)),'ko')

bar(4,mean(fit.X2Y.S2))
errorbar(4,mean(fit.X2Y.S2),std(fit.X2Y.S2)/sqrt(numel(simReps)),'ko')
bar(5,mean(fit.Y2X.S2))
errorbar(5,mean(fit.Y2X.S2),std(fit.Y2X.S2)/sqrt(numel(simReps)),'ko')

set(gca, 'XTick',[1,2,4,5], 'XTickLabel', {'S1: X->Y','S1: Y->X','S2: X->Y','S2: Y->X'})
ylabel('[bits]')
title('FIT in the two directions')
