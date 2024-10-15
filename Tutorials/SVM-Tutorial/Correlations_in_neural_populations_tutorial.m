% -------------------------------------------------------------------------
% Tutorial: Population Information Analysis
% -------------------------------------------------------------------------
% In this tutorial, we simulate and analyze population responses under
% different correlation scenarios. We perform information breakdown
% analysis of neuron pairs and SVM analysis of the whole population.

% -------------------------------------------------------------------------
% Step 1: Initialize Parameters and Options
% -------------------------------------------------------------------------

%% 5 Repetitions again!!!

clear; clc

rng(0);

nNeurons = 20; % Number of neurons
nPairs = nNeurons*(nNeurons-1)/2;
nTrials_per_stim = 400; % Number of trials per stimulus
nShuff = 2; % Number of stim-fixed shufflings per repetition
nReps = 1; % Number of repetitions of the simulation
tPoints = 1000; % ms
scenario_labels = {'noNoise','limitNoise','onlyNoise'};

% Population MI options (used for SVM)
MI_opts.verbose = false;
MI_opts.method = "dr";
MI_opts.bias = 'qe';
MI_opts.btsp = 0; % set to 0 because we implement stim-fixed shuffling later
MI_opts.bin_methodX = 'none'; % we don't need to bin the neural responses
MI_opts.bin_methodY = 'none';

% Info breakdown opts
info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
bdw_bins = 5; % number of bins used to discretize spike counts for information breakdown calculations

% Initialize structures
for scIdx = 1:numel(scenario_labels)
    scLab = scenario_labels{scIdx};
    for bdwIdx = 1:numel(info_bdw_terms)
        bdwLab = info_bdw_terms{bdwIdx};
        MI_breakdown.(scLab).(bdwLab) = nan(nReps, nPairs);
    end
    noiseCorr.(scLab) = nan(nReps, 2, nPairs);
end

% -------------------------------------------------------------------------
% Step 2: Run Simulations and Analye the Results for Each Scenario
% -------------------------------------------------------------------------

% Generate Stimulus
S = [1*ones(1,nTrials_per_stim),2*ones(1,nTrials_per_stim)];

% Here we perform target-fixed shuffling of the responses
shIdxs = nan(nShuff,nNeurons,2*nTrials_per_stim);

for shIdx = 1:nShuff
    for stimVal = 1:2
        for cellIdx = 1:nNeurons
            stim_idxs = find(S == stimVal);
            rand_idxs = randperm(numel(stim_idxs));
            val_idxs = stim_idxs(rand_idxs);
            shIdxs(shIdx,cellIdx,stim_idxs) = val_idxs;
        end
    end
end

% Loop through each scenario
for scIdx = 1:numel(scenario_labels)
    scLab = scenario_labels{scIdx};

    % ---------------------------------------------------------------------
    % Step 2.1: Simulate Spike Trains for each scenario
    % ---------------------------------------------------------------------

    switch scIdx
        case 1 % Scenario 1: uncorrelated noise
            disp(['Simulating ',scenario_labels{1},' scenario'])
            lambda_1 = 1;
            lambda_2 = 2;
            lambda_noise = 0;
            noise = false;
            [R_1, R_2, R_count_tot] = simulate_spike_trains(nNeurons, tPoints, nTrials_per_stim, lambda_1, lambda_2,lambda_noise, scLab);


        case 2 % Scenario 2: info limiting correlations
            disp(['Simulating ',scenario_labels{2},' scenario'])
            lambda_1 = 0.8;
            lambda_2 = 1.8;
            lambda_noise = 0.2;
            noise = true;
            [R_1, R_2, R_count_tot] = simulate_spike_trains(nNeurons, tPoints, nTrials_per_stim, lambda_1, lambda_2, lambda_noise, scLab);

        case 3 % Scenario 3: pure noise-corr info
            disp(['Simulating ',scenario_labels{3},' scenario'])
            lambda_1 = 1;
            lambda_2 = 2;
            lambda_noise = 1;
            noise = true;
            [R_1, R_2, R_count_tot] = simulate_spike_trains(nNeurons, tPoints, nTrials_per_stim, lambda_1, lambda_2, lambda_noise, scIdx);
    end
    

    % -------------------------------------------------------------------------
    % Step 2.2: Compute Noise Correlation and Information Breakdown
    % -------------------------------------------------------------------------

    pairIdx = 0;
    for cell1 = 1:nNeurons
        for cell2 = cell1+1:nNeurons
            pairIdx = pairIdx+1;
            jointResp = [R_count_tot(cell1,:);R_count_tot(cell2,:)];

            if ~isempty(cell2)
                noiseCorr.(scLab)(1,pairIdx) = corr(R_1(cell1,:)',R_1(cell2,:)');
                noiseCorr.(scLab)(2,pairIdx) = corr(R_2(cell1,:)',R_2(cell2,:)');
                jointResp(jointResp>bdw_bins -1)=bdw_bins -1;
                infoBdw = information(jointResp,S,MI_opts,{'I','ILIN','ISS','ICI','ICD'});
            end

            % Store information breakdown results
            for bdwIdx = 1:numel(info_bdw_terms)
                bdwLab = info_bdw_terms{bdwIdx};
                MI_breakdown.(scLab).(bdwLab)(pairIdx) = infoBdw{bdwIdx};

            end
        end
    end

    % -------------------------------------------------------------------------
    % Step 2.3: Compute Population Information with SVM Decoder
    % -------------------------------------------------------------------------

    % Set SVM Parameters
    opts.optimize_params = false;   
    opts.cv_type = "KFold"; % select method for CV partitioning (used both for information computation and for hyperparam. optimization)
    opts.K = 5; % CV fold
    opts.libsvm = false; % change to true, if you want to use libsvm instead of fitcsvm
    
    
    % Compute MI for simulated data
    % Split data into training and test
    cvPartition = cvpartition(S, opts.cv_type, opts.K); % Crossvalidation partitioning used to compute information on the held out trials
    
    % Linear SVM
    info_out_all = [];
    for i = 1:cvPartition.NumTestSets
        test_idxs = find(cvPartition.test(i));
        opts.svm_family = 'linear';
        [PredLabels] = svm_pipeline(R_count_tot',S',test_idxs, opts);
        info_out = information(PredLabels', S(test_idxs), MI_opts, {'I'});
        info_out_all(end+1) = info_out{1};
    end
    info_out_mean = mean(info_out_all);
    MI.linear.(scLab).mean = info_out_mean;
    MI.linear.(scLab).all = info_out_all;
    
    % RBF SVM
    info_out_all = [];
    for i = 1:cvPartition.NumTestSets
        test_idxs = find(cvPartition.test(i));
         opts.svm_family = 'RBF';
        [PredLabels] = svm_pipeline(R_count_tot',S',test_idxs, opts);
        info_out = information(PredLabels', S(test_idxs), MI_opts, {'I'});
        info_out_all(end+1) = info_out{1};
    end
    info_out_mean = mean(info_out_all);
    MI.RBF.(scLab).mean = info_out_mean;
    MI.RBF.(scLab).all = info_out_all;

    % Perform shufflung and compute MI for shuffled data
    for shIdx = 1:nShuff
        
        RSh_count_tot = nan(size(R_count_tot));
        for cellIdx = 1:nNeurons
            RSh_count_tot(cellIdx,:) = R_count_tot(cellIdx,squeeze(shIdxs(shIdx,cellIdx,:)));
        end

        % Linear SVM
        info_out_all = [];
        for i = 1:cvPartition.NumTestSets
            test_idxs = find(cvPartition.test(i));
            opts.svm_family = 'linear';
            [PredLabels] = svm_pipeline(RSh_count_tot',S',test_idxs, opts);
            info_out = information(PredLabels', S(test_idxs), MI_opts, {'I'});
            info_out_all(end+1) = info_out{1};
        end
        info_out_mean = mean(info_out_all);
        MISh.linear.(scLab)(shIdx).mean = info_out_mean;
        MISh.linear.(scLab)(shIdx).all = info_out_all;


        % RBF SVM
        info_out_all = [];
        for i = 1:cvPartition.NumTestSets
            test_idxs = find(cvPartition.test(i));
            opts.svm_family = 'RBF';
            [PredLabels] = svm_pipeline(RSh_count_tot',S',test_idxs, opts);
            info_out = information(PredLabels', S(test_idxs), MI_opts, {'I'});
            info_out_all(end+1) = info_out{1};
        end
        info_out_mean = mean(info_out_all);
        MISh.RBF.(scLab)(shIdx).mean = info_out_mean;
        MISh.RBF.(scLab)(shIdx).all = info_out_all;
    end
end
% Save results
save(['populationInfo_results_', date])

% -------------------------------------------------------------------------
% Step 3: Plot Results
% -------------------------------------------------------------------------
load(['populationInfo_results_', date])

maxY = 0.85;
minY = -0.05;

scenario_titles = {'no noise-corr.', 'info-limiting corr.','noise-only info.'};
figure('Position', [360, 347, 785, 271])
for scIdx = 1:numel(scenario_labels)
    scLab = scenario_labels{scIdx};
    count = 0;

    subplot(1,3,scIdx)
    hold on
    LinSh_all = [];
    RBFSh_all = [];
    for i = 1:nShuff
        LinSh_all = [LinSh_all,MISh.linear.(scLab)(i).all];
        RBFSh_all = [RBFSh_all,MISh.RBF.(scLab)(i).all];
    end 

    [meanLin, errLin] = compute_stats(MI.linear.(scLab).all);
    [meanLinSh, errLinSh] = compute_stats(LinSh_all);
    [meanRBF, errRBF] = compute_stats(MI.RBF.(scLab).all);
    [meanRBFSh, errRBFSh] = compute_stats(RBFSh_all);

    % Barplots
    h{1} = bar(scIdx,meanLin,'FaceColor',[0 0 0.7]);
    errorbar(scIdx,meanLin,errLin,'k--')
    h{2} = bar(scIdx+1,meanLinSh,'FaceColor',[0.4 0.4 0.8]);
    errorbar(scIdx+1,meanLinSh,errLinSh,'k--')
    h{3} = bar(scIdx+2,meanRBF,'FaceColor',[0.7 0 0]);
    errorbar(scIdx+2,meanRBF,errRBF,'k--')
    h{4} = bar(scIdx+3,meanRBFSh,'FaceColor',[0.8 0.4 0.4]);
    errorbar(scIdx+3,meanRBFSh,errRBFSh,'k--')

    ylabel('Information [bits]')
    if scIdx == 3
        legend([h{1} h{2} h{3} h{4}], 'Linear', 'Linear shuff.', 'RBF', 'RBF shuff.')
    end
    ylim([minY, maxY])

    title([scenario_titles{scIdx}, ' scenario'])

end
sgtitle(['Population information in different correlation scenarios, N = ', num2str(nNeurons), ' neurons'])

% -------------------------------------------------------------------------
% Step 4: Information Breakdown Plot
% -------------------------------------------------------------------------

maxY = 0.3;
minY = -0.05;
xShift = 0;
fSize = 15;
minYText = 0.03;

scenario_titles = {'no noise-corr.', 'info-limiting corr.','noise-only info.'};
figure('Position', [360, 347, 785, 271])
for scIdx = 1:numel(scenario_labels)
    scLab = scenario_labels{scIdx};
    count = 0;
    subplot(1,3,scIdx)
    hold on
   
    % Compute mean breakdown terms
    [meanJoint, errJoint] = compute_stats(MI_breakdown.(scLab).Joint);
    [meanLin, errLin] = compute_stats(MI_breakdown.(scLab).ILIN);
    [meanISS, errISS] = compute_stats(MI_breakdown.(scLab).ISS);
    [meanICI, errICI] = compute_stats(MI_breakdown.(scLab).ICI);
    [meanICD, errICD] = compute_stats(MI_breakdown.(scLab).ICD);
    
    % Barplots
    % Joint information
    h{1} = bar(scIdx,meanJoint,'FaceColor',[0.4 0.4 0.4]);
    errorbar(scIdx,meanJoint,errJoint,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).Joint,2),[],'tail','right');
    pvalues_plot(p,scIdx-xShift,max([(5/4)*meanJoint,minYText]),0.01,1,fSize,0,'k')
    % I_Lin
    h{2} = bar(scIdx+1,meanLin,'FaceColor',[0.8 0.4 0.4]);
    errorbar(scIdx+1,meanLin,errLin,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ILIN,2),[],'tail','right');
    pvalues_plot(p,scIdx+1-xShift,max([(5/4)*meanLin,minYText]),0.01,1,fSize,0,'k')
    % Sig. sim.
    h{3} = bar(scIdx+2,meanISS,'FaceColor',[0.4 0.8 0.4]);
    errorbar(scIdx+2,meanISS,errISS,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ISS,2),[],'tail','left');
    pvalues_plot(p,scIdx+2-xShift,max([(5/4)*meanISS,minYText]),0.01,1,fSize,0,'k')
    % Cor-indep
    h{4} = bar(scIdx+3,meanICI,'FaceColor',[0.4 0.4 0.8]);
    errorbar(scIdx+3,meanICI,errICI,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ICI,2));
    pvalues_plot(p,scIdx+3-xShift,max([(5/4)*meanICI,minYText]),0.01,1,fSize,0,'k')
    % Cor-dep
    h{5} = bar(scIdx+4,meanICD,'FaceColor',[0.8 0.4 0.8]);
    errorbar(scIdx+4,meanICD,errICD,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ICD,2),[],'tail','right');
    pvalues_plot(p,scIdx+4-xShift,max([(5/4)*meanICD,minYText]),0.01,1,fSize,0,'k')

    ylabel('Information [bits]')
    if scIdx == 3
        legend([h{1} h{2} h{3} h{4} h{5}], 'Joint', 'Linear', 'Sig. sim.', 'Corr. indep.', 'Corr. dep.')
    end
    ylim([minY,maxY])

    title([scenario_titles{scIdx}, ' scenario'])
end
sgtitle(['Information breakdown in different correlation scenarios, N = ',num2str(nNeurons),' neurons'])


% -------------------------------------------------------------------------
% Helper Function: Simulate Spike Trains
% -------------------------------------------------------------------------
function [mean_val, sem_val] = compute_stats(data)
    mean_val = mean(data);
    std_val = std(data);
    n = length(data);
    sem_val = std_val / sqrt(n);
end

function [R_1, R_2, R_count_tot] = simulate_spike_trains(nNeurons, tPoints, nTrials_per_stim, lambda_1, lambda_2, lambda_noise, scIdx)
R_1 = zeros(nNeurons,tPoints,nTrials_per_stim);
R_2 = zeros(nNeurons,tPoints,nTrials_per_stim);

if strcmp(scIdx, 'noNoise')
    for cellIdx = 1:nNeurons
        for trialIdx = 1:nTrials_per_stim
            [spikes1] = poisson_spike_gen(1:tPoints, lambda_1/tPoints, 0);
            [spikes2] = poisson_spike_gen(1:tPoints, lambda_2/tPoints, 0);
            R_1(cellIdx,:,trialIdx) = spikes1;
            R_2(cellIdx,:,trialIdx) = spikes2;
        end
    end
else
    for trialIdx = 1:nTrials_per_stim
        [spikes1_noise] = poisson_spike_gen(1:tPoints, lambda_noise/tPoints, 0);
        if strcmp(scIdx, 'limitNoise')
            [spikes2_noise] = poisson_spike_gen(1:tPoints, lambda_noise/tPoints, 0);
        end
        for cellIdx = 1:nNeurons
            [spikes1] = poisson_spike_gen(1:tPoints, lambda_1/tPoints, 0);
            [spikes2] = poisson_spike_gen(1:tPoints, lambda_2/tPoints, 0);
            R_1(cellIdx,:,trialIdx) = spikes1 + spikes1_noise;
            if strcmp(scIdx, 'limitNoise')
                R_2(cellIdx,:,trialIdx) = spikes2 + spikes2_noise;
            else
                R_2(cellIdx,:,trialIdx) = spikes2;
            end
        end
    end
end
% Compute single-trial features
R_1_count = squeeze(sum(R_1,2));
R_2_count = squeeze(sum(R_2,2));
R_count_tot = [R_1_count,R_2_count];

end


