% This is a tutorial demonstrating the use of functions for calculating 
% FIT, and transfer entropy TE in simulated EEG signals across multiple populations and trials. 
% The simulation generates neural signals, computes mutual information, 
% and evaluates the effects of stimuli while performing constrained and 
% unconstrained shuffling for null hypothesis testing.

clear all; clc;
rng(0); 

% Simulation parameters
simReps = 5; % Number of simulation repetitions
nPopulations = 4; % Number of populations (e.g., brain regions)
nSubpopulations = 2; % Number of subpopulations within each population
nTrials = 200; % Number of trials
nTimepoints = 30; % Number of time points in each trial
delay = 5; % Delay in time points for information transfer between populations
stim_window = [3 12]; % Time window for the stimuli input
nShuff = 10; % Number of shuffles for null hypothesis testing
null_samples = 500; % Number of null samples to generate

c = 1; % Coefficient for coupling between populations

% Initialize variables for storing signal data and results
X = zeros(nPopulations, nSubpopulations, nTimepoints, nTrials);
infoS1 = zeros(nPopulations, nTimepoints, simReps);
infoS2 = zeros(nPopulations, nTimepoints, simReps);

infoS1_sh = zeros(nPopulations, nTimepoints, simReps, nShuff);
infoS2_sh = zeros(nPopulations, nTimepoints, simReps, nShuff);

% Initialize matrices for FIT and TE calculations
FIT1 = zeros(nPopulations, nPopulations, simReps);
FIT2 = zeros(nPopulations, nPopulations, simReps);
cFIT1 = zeros(nPopulations, nPopulations, simReps);
cFIT2 = zeros(nPopulations, nPopulations, simReps);
TE = zeros(nPopulations, nPopulations, simReps);

% Initialize options for information calculation
opts.verbose = false; % Set verbosity of output
opts.method = "dr"; % Method for estimating information
opts.bias = 'naive'; % Bias correction method
opts.xtrp = 10; % Extrapolation parameter
opts.btsp = 0; % Bootstrap option (0 = off)
opts.n_binsS = 2; % Number of Stimuli
opts.n_binsX = 3; % Number of bins for X data
opts.n_binsY = 3; % Number of bins for Y data
opts.bin_methodS = 'none'; % No binning method for stimulus
opts.bin_methodX = 'eqpop'; % Equal populated binning for X data
opts.taux = -delay; % Time lag for X
opts.tauy = -delay; % Time lag for Y
opts.tpres = 7 + delay;
opts.parallel = 0; % Parallel processing option (0 = off)
opts.singleTimepoint = true; % Single time point option for TE anaylsis

% Simulate signal for each subpopulation over time
for repIdx = 1:simReps
    disp(repIdx); 
    X = zeros(nPopulations, nSubpopulations, nTimepoints, nTrials); % Reset X for each repetition
    S1 = randi(opts.n_binsS, 1, nTrials); % Random stimulus assignments for S1
    S2 = randi(opts.n_binsS, 1, nTrials); % Random stimulus assignments for S2

    epsX = 0.5; % Noise level for signals

    stim11 = stim_window(1); % Start of stimulus input
    stim12 = stim_window(2); % End of stimulus input
    
    % Generate neural signal over time
    for t = 1:nTimepoints
        X(:,:,t,:) = epsX * randn(nPopulations, nSubpopulations, 1, nTrials); % Initialize noise

        % Add stimulus effects
        if (t > stim11) && (t < stim12)
            X(1, 1, t, :) = squeeze(X(1, 1, t, :)) + epsX * randn + 2 * (S1 - 1.5)'; % Modify population 1 for S1
            X(2, 2, t, :) = squeeze(X(2, 2, t, :)) + epsX * randn + 2 * (S2 - 1.5)'; % Modify population 2 for S2
        end

        % Information transfer between the subpopulations
        if t > delay
            X(1, 2, t, :) = X(1, 2, t, :) + c * X(2, 2, t - delay, :); 
            X(2, 1, t, :) = X(2, 1, t, :) + c * X(1, 1, t - delay, :); 
            X(3, 1, t, :) = X(3, 1, t, :) + c * X(1, 1, t - delay, :); 
            X(4, 2, t, :) = X(4, 1, t, :) + c * X(3, 1, t - delay, :); 
        end

        % Calculate mutual information (MI) for both stimuli
        opts.bin_methodY = 'none'; % No binning method for Y
        for neuron = 1:nPopulations
            infoS1(neuron, t, repIdx) = cell2mat(information(squeeze(sum(X(neuron, :, t, :), 2))', S1, opts, {'I'})); % MI for S1
            infoS2(neuron, t, repIdx) = cell2mat(information(squeeze(sum(X(neuron, :, t, :), 2))', S2, opts, {'I'})); % MI for S2
        end
    end
    
    X = squeeze(sum(X, 2)); % Collapse the second dimension

    % Compute temporal profile of FIT for each pair of populations
    opts.bin_methodY = 'eqpop'; 
    t = 7 + delay; % Current time point for analysis, we choose 7 as it is after the stimulus input 
    
    % Compute FIT and TE for all pairs of populations
    for pop1 = 1:nPopulations
        for pop2 = 1:nPopulations
            if pop1 ~= pop2
              
                FIT1(pop1, pop2, repIdx) = FIT(S1, X(pop1,:, :), X(pop2,:, :), opts); % FIT from pop1 to pop2 with S1              
                FIT2(pop1, pop2, repIdx) = FIT(S2, X(pop1,:, :), X(pop2,:, :), opts); % FIT from pop2 to pop1 with S2
                    %cFIT(pop1, pop2, repIdx) = cFIT()
                TE(pop1, pop2, repIdx) = cell2mat(transferentropy2(X(pop1, 1:t, :), X(pop2, 1:t, :), opts, {'TE'})); % Transfer entropy from pop1 to pop2
            end
        end
    end

    % Compute null hypothesis for MI, FIT, and TE through shuffling
    for shIdx = 1:nShuff
        disp(shIdx); % Display the current shuffle index
        shuffledS1 = S1(randperm(numel(S1))); % Shuffle S1
        shuffledS2 = S2(randperm(numel(S2))); % Shuffle S2
        S1_S2 = map_Nd_array_to_1d([S1; S2]); % Map S1 and S2 into a 1D array
        for pop1 = 1:nPopulations
            % Permute data to have trials in the first dimension
            sender = X(pop1, :, :); % Extract sender population data
            sender_perm = permute(sender, [3, 1, 2]); % Permute dimensions for shuffling

            % Constrained shuffling of fixed values of S1 and S2
            sender_s1_s2_shuff_perm = shuffle(S1_S2', sender_perm, 1, [1, 0, 0]); % Constrained shuffle
            sender_s1_s2_shuff = permute(sender_s1_s2_shuff_perm, [2, 3, 1]); % Restore original dimensions

            % Unconstrained shuffling of the sender
            sender_all_shuff_perm = shuffle(S1', sender_perm, 0, [1, 0, 0]); % Unconstrained shuffle
            sender_all_shuff = permute(sender_all_shuff_perm, [2, 3, 1]); % Restore original dimensions

            % Shuffle the sender at fixed values of S1
            sender_s1_shuff_perm = shuffle(S1', sender_perm, 1, [1, 0, 0]); % Fixed S1 shuffle
            sender_s1_shuff = permute(sender_s1_shuff_perm, [2, 3, 1]); % Restore original dimensions

            % Shuffle the sender at fixed values of S2
            sender_s2_shuff_perm = shuffle(S2', sender_perm, 1, [1, 0, 0]); % Fixed S2 shuffle
            sender_s2_shuff = permute(sender_s2_shuff_perm, [2, 3, 1]); % Restore original dimensions

            % Initialize shuffled data
            X_s1_shuff = X; X_s2_shuff = X; X_all_shuff = X; X_s1_s2_shuff = X;
            X_s1_shuff(pop1, :, :) = sender_s1_shuff; % Replace with shuffled S1 data
            X_s2_shuff(pop1, :, :) = sender_s2_shuff; % Replace with shuffled S2 data
            X_s1_s2_shuff(pop1, :, :) = sender_s1_s2_shuff; % Replace with shuffled S1/S2 data
            X_all_shuff(pop1, :, :) = sender_all_shuff; % Replace with shuffled all data

            t = 7 + delay; % Current time point for analysis
            X_pres = squeeze(X(:, t, :)); % Current state of populations
            X_past = squeeze(X(:, t - delay, :)); % Previous state of populations

            % Compute past and present of neural signals for shuffled data
            X_s1_shuff_pres = squeeze(X_s1_shuff(:, t, :)); % Present state of shuffled S1
            X_s1_shuff_past = squeeze(X_s1_shuff(:, t - delay, :)); % Past state of shuffled S1
            X_s2_shuff_pres = squeeze(X_s2_shuff(:, t, :)); % Present state of shuffled S2
            X_s2_shuff_past = squeeze(X_s2_shuff(:, t - delay, :)); % Past state of shuffled S2

            for pop2 = 1:nPopulations
                if pop1 ~= pop2
                    % Compute FIT and TE for shuffled data
                    FIT1_sh.conditioned(pop1, pop2, repIdx, shIdx) = FIT(S1, X_s1_shuff_past(pop1, :), X_s1_shuff_past(pop2, :), X_s1_shuff_pres(pop2, :), opts); % FIT conditioned on S1
                    FIT2_sh.conditioned(pop1, pop2, repIdx, shIdx) = FIT(S2, X_s2_shuff_past(pop1, :), X_s2_shuff_past(pop2, :), X_s2_shuff_pres(pop2, :), opts); % FIT conditioned on S2

                    % Simple S shuffling
                    FIT1_sh.simple(pop1, pop2, repIdx, shIdx) = FIT(shuffledS1, X_past(pop1, :), X_past(pop2, :), X_pres(pop2, :), opts); % FIT with shuffled S1
                    FIT2_sh.simple(pop1, pop2, repIdx, shIdx) = FIT(shuffledS2, X_past(pop1, :), X_past(pop2, :), X_pres(pop2, :), opts); % FIT with shuffled S2

                    % Compute transfer entropy for shuffled data
                    TE_sh.simple(pop1, pop2, repIdx, shIdx) = cell2mat(transferentropy2(X_all_shuff(pop1, 1:t, :), X_all_shuff(pop2, 1:t, :), opts, {'TE'})); % TE for all shuffled data
                    TE_sh.conditioned(pop1, pop2, repIdx, shIdx) = cell2mat(transferentropy2(X_s1_s2_shuff(pop1, 1:t, :), X_s1_s2_shuff(pop2, 1:t, :), opts, {'TE'})); % TE for conditioned shuffled data
                end
            end
        end
    end
end
