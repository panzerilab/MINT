% Seed for reproducibility
rng('default');

% Time and trial settings
dt = 1/500;
trialendtime = 0.4;
t_trial = 0:dt:trialendtime;
nStimuli = 2;

% Number of repetitions
nRepetitions = 4;

% Mutual Information options
shuffopts.method = 'dr';
shuffopts.bias = 'naive';
shuffopts.bin_methodX = 'eqpop';
shuffopts.bin_methodY = 'none';
shuffopts.n_binsX = 2;
shuffopts.n_binsY = 0;
shuffopts.btsp = 10;
shuffopts.compute_nulldist = false;
shuffopts.parallel = 0;
shuffopts.supressWarnings = true;
shuffopts.verbose = false;

qeopts = shuffopts;
qeopts.bias = 'qe';
qeopts.xtrp = 10;
qeopts.btsp = 0;

shopts = shuffopts;
shopts.bias = 'naive';
shopts.btsp = 0;

ptopts = shuffopts;
ptopts.bias = 'pt';
ptopts.xtrp = 0;
ptopts.btsp = 0;

bubopts = shuffopts;
bubopts.bias = 'bub';
bubopts.xtrp = 0;
bubopts.btsp = 0;

newshuffopts.bias = 'shuffSub';
newshuffopts.bin_method = {'eqpop', 'none'};
newshuffopts.n_bins = {2, 0};
newshuffopts.shuff = 10;
newshuffopts.compute_nulldist = false;
newshuffopts.parallel = 0;
newshuffopts.supressWarnings = true;

newqeopts = shuffopts;
newqeopts.bias = 'qe';
newqeopts.xtrp = 10;
newqeopts.shuff = 0;

newshopts = shuffopts;
newshopts.bias = 'shuffCorr';
newshopts.xtrp = 0;
newshopts.shuff = 0;

newptopts = shuffopts;
newptopts.bias = 'pt';
newptopts.xtrp = 0;
newptopts.shuff = 0;

newbubopts = shuffopts;
newbubopts.bias = 'bub';
newbubopts.xtrp = 0;
newbubopts.shuff = 0;


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
MI_naive = zeros(2,nRepetitions, length(nTrials_list));
MI_shuff = zeros(2,nRepetitions, length(nTrials_list));
MI_qe    = zeros(2,nRepetitions, length(nTrials_list));
MI_sh    = zeros(2,nRepetitions, length(nTrials_list));
MI_pt    = zeros(2,nRepetitions, length(nTrials_list));
MI_bub   = zeros(2,nRepetitions, length(nTrials_list));

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
        [mi_naive]           = information(R,S, shuffopts, {'I'});
        mi_qe                = information(R,S, qeopts, {'I'});
        mi_sh                = information(R,S, shopts, {'Ish'});
        mi_pt                = information(R,S, ptopts, {'I'});
        mi_bub               = information(R,S, bubopts, {'I'});

        [newmi_shuff, newmi_naive] = MI({R;S}, {'I(A;B)'}, newshuffopts);
        newmi_qe                   = MI({R;S}, {'I(A;B)'}, newqeopts);
        newmi_sh                   = MI({R;S}, {'I(A;B)'}, newshopts);
        newmi_pt                   = MI({R;S}, {'I(A;B)'}, newptopts);
        newmi_bub                  = MI({R;S}, {'I(A;B)'}, newbubopts);

        % Store results
        MI_naive(:, rep, tidx) = [mi_naive{1}(1), newmi_naive{1}];
        MI_shuff(:,rep, tidx)  = [mi_naive{1}(1) - mean(mi_naive{1}(2:end)), newmi_shuff{1}];
        MI_qe(:,rep, tidx)     = [mi_qe{1}, newmi_qe{1}];
        MI_sh(:,rep, tidx)     = [mi_sh{1}, newmi_sh{1}];
        MI_pt(:,rep, tidx)     = [mi_pt{1}, newmi_pt{1}];
        MI_bub(:,rep, tidx)    = [mi_bub{1}, newmi_bub{1}];
    end

end



%%

% Average MI values across 50 repetitions
MI_naive_avg = mean(MI_naive, 2);
MI_shuff_avg = mean(MI_shuff, 2);
MI_qe_avg    = mean(MI_qe, 2);
MI_sh_avg    = mean(MI_sh, 2);
MI_pt_avg    = mean(MI_pt, 2);
MI_bub_avg   = mean(MI_bub, 2);

%%
% Plotting
figure;
plot(5:9, MI_naive_avg(1,:), '-o', 'LineWidth', 2);
hold on;
plot(5:9, MI_shuff_avg(1,:), '-o', 'LineWidth', 2);
plot(5:9, MI_qe_avg(1,:),    '-o', 'LineWidth', 2);
plot(5:9, MI_sh_avg(1,:),    '-o', 'LineWidth', 2);
plot(5:9, MI_pt_avg(1,:),    '-o', 'LineWidth', 2);
plot(5:9, MI_bub_avg(1,:),   '-o', 'LineWidth', 2);

legend({'MI_{naive}', 'MI_{shuff}', 'MI_{qe}', 'MI_{sh}', 'MI_{pt}', 'MI_{bub}'}, 'Location', 'best');
xlabel('Sample Size (Number of Trials)');
ylabel('Average Mutual Information (MI)');
title('Average old MI Across Sample Sizes');

figure;
plot(5:9, MI_naive_avg(2,:), '--', 'LineWidth', 2);
hold on;
plot(5:9, MI_shuff_avg(2,:), '--', 'LineWidth', 2);
plot(5:9, MI_qe_avg(2,:),    '--', 'LineWidth', 2);
plot(5:9, MI_sh_avg(2,:),    '--', 'LineWidth', 2);
plot(5:9, MI_pt_avg(2,:),    '--', 'LineWidth', 2);
plot(5:9, MI_bub_avg(2,:),   '--', 'LineWidth', 2);

legend({'MI_{naive}', 'MI_{shuff}', 'MI_{qe}', 'MI_{sh}', 'MI_{pt}', 'MI_{bub}'}, 'Location', 'best');
xlabel('Sample Size (Number of Trials)');
ylabel('Average Mutual Information (MI)');
title('Average new MI Across Sample Sizes');
