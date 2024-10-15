% Seed for reproducibility
rng('default');

% Time and trial settings
dt = 1/500;
trialendtime = 0.4;
t_trial = 0:dt:trialendtime;
nStimuli = 2;

% Number of repetitions
nRepetitions = 10;

% Mutual Information options
shuffopts.bias = 'shuffSub';
shuffopts.bin_method = {'eqpop', 'none'};
shuffopts.n_bins = {2, 0};
shuffopts.shuff = 10;
shuffopts.compute_nulldist = false;
shuffopts.parallel = 0;
shuffopts.supressWarnings = true;

qeopts = shuffopts;
qeopts.bias = 'qe';
qeopts.xtrp = 20;
qeopts.shuff = 0;

shopts = shuffopts;
shopts.bias = 'shuffCorr';
shopts.xtrp = 0;
shopts.shuff = 0;

ptopts = shuffopts;
ptopts.bias = 'pt';
ptopts.xtrp = 0;
ptopts.shuff = 0;

bubopts = shuffopts;
bubopts.bias = 'bub';
bubopts.xtrp = 0;
bubopts.shuff = 0;

% List of trial counts to iterate over
nTrials_list = 2.^ [5:9]; %[8, 16, 32, 64, 256, 512, 1024];

% Generate response to stimuli
maxrate = [50 20];
rate = nan(nStimuli, length(t_trial));
for i = 1:nStimuli
    signal = normpdf(t_trial, 0.2, 0.05);
    rate(i,:) = maxrate(i) * signal / max(signal);
end

% Initialize matrices to store MI results for each repetition
MI_naive = zeros(nRepetitions, length(nTrials_list));
MI_shuff = zeros(nRepetitions, length(nTrials_list));
MI_qe    = zeros(nRepetitions, length(nTrials_list));
MI_sh    = zeros(nRepetitions, length(nTrials_list));
MI_pt    = zeros(nRepetitions, length(nTrials_list));
MI_bub   = zeros(nRepetitions, length(nTrials_list));

% Loop over sample sizes
for tidx = 1:length(nTrials_list)
    nTrials = nTrials_list(tidx);
    sprintf('Calculating sample size %d', nTrials)

    parfor rep = 1:nRepetitions
        % Generate responses (spike counts)
        R = [];
        S = [];
        for i = 1:nStimuli
            rate_stim = rate(i,:);
            for j = 1:nTrials
                spike_train1 = poisson_spike_gen(t_trial, rate_stim, 0);
                spike_train2 = poisson_spike_gen(t_trial, rate_stim, 0);
                % spike_train3 = poisson_spike_gen(t_trial, rate_stim, 0);
                r_tot = [sum(spike_train1); sum(spike_train2)];
                R = [R, r_tot];
                S = [S i];
            end
        end

        % Calculate Mutual Information using different methods
        [mi_shuff, mi_naive] = MI({R;S}, {'I(A;B)'}, shuffopts);
        mi_qe                = MI({R;S}, {'I(A;B)'}, qeopts);
        mi_sh                = MI({R;S}, {'I(A;B)'}, shopts);
        mi_pt                = MI({R;S}, {'I(A;B)'}, ptopts);
        mi_bub               = MI({R;S}, {'I(A;B)'}, bubopts);

        % Store results
        MI_naive(rep, tidx) = mi_naive{1};
        MI_shuff(rep, tidx) = mi_shuff{1};
        MI_qe(rep, tidx)    = mi_qe{1};
        MI_sh(rep, tidx)    = mi_sh{1};
        MI_pt(rep, tidx)    = mi_pt{1};
        MI_bub(rep, tidx)   = mi_bub{1};
    end

end

% Average MI values across 50 repetitions
MI_naive_avg = mean(MI_naive, 1);
MI_shuff_avg = mean(MI_shuff, 1);
MI_qe_avg    = mean(MI_qe, 1);
MI_sh_avg    = mean(MI_sh, 1);
MI_pt_avg    = mean(MI_pt, 1);
MI_bub_avg   = mean(MI_bub, 1);

%%
% Plotting
figure;
semilogx(nTrials_list, MI_naive_avg, '-o', 'LineWidth', 2);
hold on;
semilogx(nTrials_list, MI_shuff_avg, '-o', 'LineWidth', 2);
semilogx(nTrials_list, MI_qe_avg, '-o', 'LineWidth', 2);
semilogx(nTrials_list, MI_sh_avg, '-o', 'LineWidth', 2);
semilogx(nTrials_list, MI_pt_avg, '-o', 'LineWidth', 2);
semilogx(nTrials_list, MI_bub_avg, '-o', 'LineWidth', 2);

legend({'MI_{naive}', 'MI_{shuff}', 'MI_{qe}', 'MI_{sh}', 'MI_{pt}', 'MI_{bub}'}, 'Location', 'best');
xlabel('Sample Size (Number of Trials)');
ylabel('Average Mutual Information (MI)');
title('Average MI Across Sample Sizes');
grid on;
