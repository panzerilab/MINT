% Population information simulation and analysis
% Here we combine information breakdown analysis of neurons pairs and
% SVM analysis of the whole population

clear all; clc

rng(0);

nNeurons = 20; % Number of neurons
nPairs = nNeurons*(nNeurons-1)/2;
nTrials_per_stim = 400; % Number of trials per stimulus
nShuff = 1; % Number of stim-fixed shufflings per repetition (used to statistically destroy noise correlations)
nReps = 5; % Number of repetitions of the simulation
tPoints = 1000; % ms 
scenario_labels = {'noNoise','limitNoise','onlyNoise'};

% Population MI options (used for SVM)
MI_opts.verbose = false;
MI_opts.method = "dr";
MI_opts.bias = 'qe';
MI_opts.shuff = 0; % set to 0 because we implement stim-fixed shuffling later in the script
MI_opts.bin_methodX = 'none'; % we don't need to bin the neural responses, nor the stimulus, for the SVM analysis
MI_opts.bin_methodY = 'none';

% Info breakdown opts
info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
bdw_bins = 5; % number of binf used to discretize spike counts for information breakdown calculations

% SVM options
cv_type = "KFold"; % select method for CV partitioning (used both for information computation and for hyperparam. optimization)
folds = 5; % CV folds
optimize_params = false; % here we don't optimize hyperparameters to speed up the computation (we use the default libsvm ones)

% Initialize structures
for scIdx = 1:numel(scenario_labels)
    scLab = scenario_labels{scIdx};
    for bdwIdx = 1:numel(info_bdw_terms)
        bdwLab = info_bdw_terms{bdwIdx};
        MI_breakdown.(scLab).(bdwLab) = nan(nReps, nPairs); 
    end
    noiseCorr.(scLab) = nan(nReps, 2, nPairs);
end

%% Run nReps repetitions of the three scenarios

for repIdx = 1:nReps
    disp(['Repetition number ', num2str(repIdx)])
    S = [1*ones(1,nTrials_per_stim),2*ones(1,nTrials_per_stim)];

    % Here we perform target-fixed shuffling of the responses to destroy 
    % noise correlations and have an indep-cells null hypothesis
    shIdxs = nan(nShuff,nNeurons,2*nTrials_per_stim); % dimension 2 is the cell --> we have to permute the two neural responses independently while keeping the marginals of individual cells with target fixed

    for shIdx = 1:nShuff
        for stimVal = 1:2
            for cellIdx = 1:nNeurons
                stim_idxs = find(S == stimVal);
                rand_idxs = randperm(numel(stim_idxs));
                val_idxs = stim_idxs(rand_idxs); % permute trials at fixed target value
                shIdxs(shIdx,cellIdx,stim_idxs) = val_idxs;
            end
        end
    end

    %% Scenario 1: uncorrelated noise
    lambda_1 = 1;
    lambda_2 = 2;

    disp(['Simulating ',scenario_labels{1},' scenario'])
    % Simulate spike trains and calcium traces
    R_1 = zeros(nNeurons,tPoints,nTrials_per_stim); 
    R_2 = zeros(nNeurons,tPoints,nTrials_per_stim); 
    for cellIdx = 1:nNeurons
        for trialIdx = 1:nTrials_per_stim
            [spikes1] = poisson_spike_gen(1:tPoints, lambda_1/tPoints, 0);
            [spikes2] = poisson_spike_gen(1:tPoints, lambda_2/tPoints, 0);
            R_1(cellIdx,:,trialIdx) = spikes1;
            R_2(cellIdx,:,trialIdx) = spikes2;
        end
    end
    
    % Compute single-trial features
    R_1_count = squeeze(sum(R_1,2));
    R_2_count = squeeze(sum(R_2,2));
    R_count_tot = [R_1_count,R_2_count];
    
    % Compute noise corr and information breakdown
    pairIdx = 0;
    for cell1 = 1:nNeurons
        for cell2 = cell1+1:nNeurons
            pairIdx = pairIdx+1;
            jointResp = [R_count_tot(cell1,:);R_count_tot(cell2,:)];
            
            if ~isempty(cell2)
                noiseCorr.noNoise(repIdx,1,pairIdx) = corr(R_1(cell1,:)',R_1(cell2,:)');
                noiseCorr.noNoise(repIdx,2,pairIdx) = corr(R_2(cell1,:)',R_2(cell2,:)');
                
                % Set responses with more than bdw_bins-1 spikes to bdw_bins-1
                jointResp(jointResp>bdw_bins -1)=bdw_bins -1;

                % We use MI_opts for electrophys recordings since we don't
                % want to bin
                infoBdw = information(jointResp,S,MI_opts,{'I','ILIN','ISS','ICI','ICD'});
            end
            
            MI_breakdown.noNoise.Joint(repIdx,pairIdx) = infoBdw{1}; 
            MI_breakdown.noNoise.ILIN(repIdx,pairIdx) = infoBdw{2}; 
            MI_breakdown.noNoise.ISS(repIdx,pairIdx) = infoBdw{3}; 
            MI_breakdown.noNoise.ICI(repIdx,pairIdx) = infoBdw{4}; 
            MI_breakdown.noNoise.ICD(repIdx,pairIdx) = infoBdw{5}; 
            
        end
    end

    % Compute population info
    [MI.linear.noNoise(repIdx)] = svm_MI(MI_opts, R_count_tot', S', cv_type, 'linear_SVM', folds, optimize_params);
    [MI.RBF.noNoise(repIdx)] = svm_MI(MI_opts, R_count_tot', S', cv_type, 'RBF_SVM', folds, optimize_params);
    [MI.pca.noNoise(repIdx)] = pca_MI(MI_opts, R_count_tot, S);
    [MI.nnmf.noNoise(repIdx)] = nnmf_MI(MI_opts, R_count_tot, S);
    
    for shIdx = 1:nShuff
        RSh_count_tot = nan(size(R_count_tot));

        for cellIdx = 1:nNeurons
            RSh_count_tot(cellIdx,:) = R_count_tot(cellIdx,squeeze(shIdxs(shIdx,cellIdx,:)));
        end
        [MISh.linear.noNoise(repIdx,shIdx)] = svm_MI(MI_opts, RSh_count_tot', S', cv_type, 'linear_SVM', folds, optimize_params);
        [MISh.RBF.noNoise(repIdx,shIdx)] = svm_MI(MI_opts, RSh_count_tot', S', cv_type, 'RBF_SVM', folds, optimize_params);
        [MISh.pca.noNoise(repIdx,shIdx)] = pca_MI(MI_opts, RSh_count_tot, S);
        [MISh.nnmf.noNoise(repIdx,shIdx)] = nnmf_MI(MI_opts, RSh_count_tot, S);
    end
    %% Scenario 2: info limiting correlations
    disp(['Simulating ',scenario_labels{2},' scenario'])
    
    lambda_1 = 0.8;
    lambda_2 = 1.8;
    lambda_noise = 0.2;

    % Simulate spike trains and calcium traces
    R_1 = zeros(nNeurons,tPoints,nTrials_per_stim); 
    R_2 = zeros(nNeurons,tPoints,nTrials_per_stim);

    for trialIdx = 1:nTrials_per_stim
        % shared noise: we generate it separately for stim1 and stim2 simply 
        % because we are looping within trials per stimulus, nonetheless the
        % properties of the shared noise are the same for the two stimuli
        [spikes1_noise] = poisson_spike_gen(1:tPoints, lambda_noise/tPoints, 0);
        [spikes2_noise] = poisson_spike_gen(1:tPoints, lambda_noise/tPoints, 0);

        for cellIdx = 1:nNeurons
            [spikes1] = poisson_spike_gen(1:tPoints, lambda_1/tPoints, 0);
            [spikes2] = poisson_spike_gen(1:tPoints, lambda_2/tPoints, 0);
            R_1(cellIdx,:,trialIdx) = spikes1 + spikes1_noise;
            R_2(cellIdx,:,trialIdx) = spikes2 + spikes2_noise;
        end
    end

    % Compute single-trial features
    R_1_count = squeeze(sum(R_1,2));
    R_2_count = squeeze(sum(R_2,2));

    R_count_tot = [R_1_count,R_2_count];
    
    % Compute noise corr and information breakdown
    pairIdx = 0;
    for cell1 = 1:nNeurons
        for cell2 = cell1+1:nNeurons
            pairIdx = pairIdx+1;
            jointResp = [R_count_tot(cell1,:);R_count_tot(cell2,:)];
            
            if ~isempty(cell2)
                noiseCorr.limitNoise(repIdx,1,pairIdx) = corr(R_1(cell1,:)',R_1(cell2,:)');
                noiseCorr.limitNoise(repIdx,2,pairIdx) = corr(R_2(cell1,:)',R_2(cell2,:)');
                
                jointResp(jointResp>bdw_bins -1)=bdw_bins -1;

                infoBdw = information(jointResp,S,MI_opts,{'I','ILIN','ISS','ICI','ICD'});
            end
            
            MI_breakdown.limitNoise.Joint(repIdx,pairIdx) = infoBdw{1}; 
            MI_breakdown.limitNoise.ILIN(repIdx,pairIdx) = infoBdw{2}; 
            MI_breakdown.limitNoise.ISS(repIdx,pairIdx) = infoBdw{3}; 
            MI_breakdown.limitNoise.ICI(repIdx,pairIdx) = infoBdw{4}; 
            MI_breakdown.limitNoise.ICD(repIdx,pairIdx) = infoBdw{5};

        end
    end
    
    % Compute population info
    [MI.linear.limitNoise(repIdx)] = svm_MI(MI_opts, R_count_tot', S', cv_type, 'linear_SVM', folds, optimize_params);
    [MI.RBF.limitNoise(repIdx)] = svm_MI(MI_opts, R_count_tot', S', cv_type, 'RBF_SVM', folds, optimize_params);
    [MI.pca.limitNoise(repIdx)] = pca_MI(MI_opts, R_count_tot, S);
    [MI.nnmf.limitNoise(repIdx)] = nnmf_MI(MI_opts, R_count_tot, S);
    for shIdx = 1:nShuff 
        RSh_count_tot = nan(size(R_count_tot));
        for cellIdx = 1:nNeurons
            RSh_count_tot(cellIdx,:) = R_count_tot(cellIdx,squeeze(shIdxs(shIdx,cellIdx,:)));
        end
        [MISh.linear.limitNoise(repIdx,shIdx)] = svm_MI(MI_opts, RSh_count_tot', S', cv_type, 'linear_SVM', folds, optimize_params);
        [MISh.RBF.limitNoise(repIdx,shIdx)] = svm_MI(MI_opts, RSh_count_tot', S', cv_type, 'RBF_SVM', folds, optimize_params);
        [MI.pca.limitNoise(repIdx,shIdx)] = pca_MI(MI_opts, RSh_count_tot, S);
        [MI.nnmf.limitNoise(repIdx,shIdx)] = nnmf_MI(MI_opts, RSh_count_tot, S);
    end
    %% Scenario 3: pure noise-corr info
    disp(['Simulating ',scenario_labels{3},' scenario'])
    
    lambda_1 = 1;
    lambda_2 = 2;
    lambda_noise1 = 1;

    % Simulate spike trains and calcium traces
    R_1 = zeros(nNeurons,tPoints,nTrials_per_stim); R_1_Ca = R_1;
    R_2 = zeros(nNeurons,tPoints,nTrials_per_stim); R_2_Ca = R_2;

    for trialIdx = 1:nTrials_per_stim
        % shared noise: we generate it separately for stim1 and stim2 simply 
        % because we are looping within trials per stimulus, nonetheless the
        % properties of the shared noise are the same for the two stimuli
        [spikes1_noise] = poisson_spike_gen(1:tPoints, lambda_noise1/tPoints, 0);

        for cellIdx = 1:nNeurons
            [spikes1] = poisson_spike_gen(1:tPoints, lambda_1/tPoints, 0);
            [spikes2] = poisson_spike_gen(1:tPoints, lambda_2/tPoints, 0);
            R_1(cellIdx,:,trialIdx) = spikes1 + spikes1_noise;
            R_2(cellIdx,:,trialIdx) = spikes2;
        end
    end

    % Compute single-trial features
    R_1_count = squeeze(sum(R_1,2));
    R_2_count = squeeze(sum(R_2,2));

    R_count_tot = [R_1_count,R_2_count];

    % Compute noise corr and information breakdown
    pairIdx = 0;
    for cell1 = 1:nNeurons
        for cell2 = cell1+1:nNeurons
            pairIdx = pairIdx+1;
            jointResp = [R_count_tot(cell1,:);R_count_tot(cell2,:)];

            if ~isempty(cell2)
                noiseCorr.onlyNoise(repIdx,1,pairIdx) = corr(R_1(cell1,:)',R_1(cell2,:)');
                noiseCorr.onlyNoise(repIdx,2,pairIdx) = corr(R_2(cell1,:)',R_2(cell2,:)');
                
                jointResp(jointResp>bdw_bins -1)=bdw_bins -1;

                % We use MI_opts for electrophys recordings since we don't
                % want to bin
                infoBdw = information(jointResp,S,MI_opts,{'I','ILIN','ISS','ICI','ICD'});
            end
            
            MI_breakdown.onlyNoise.Joint(repIdx,pairIdx) = infoBdw{1}; 
            MI_breakdown.onlyNoise.ILIN(repIdx,pairIdx) = infoBdw{2};
            MI_breakdown.onlyNoise.ISS(repIdx,pairIdx) = infoBdw{3}; 
            MI_breakdown.onlyNoise.ICI(repIdx,pairIdx) = infoBdw{4}; 
            MI_breakdown.onlyNoise.ICD(repIdx,pairIdx) = infoBdw{5}; 
            
        end
    end
    
    % Compute population info
    [MI.linear.onlyNoise(repIdx)] = svm_MI(MI_opts, R_count_tot', S', cv_type, 'linear_SVM', folds, optimize_params);
    [MI.RBF.onlyNoise(repIdx)] = svm_MI(MI_opts, R_count_tot', S', cv_type, 'RBF_SVM', folds, optimize_params);
    [MI.pca.onlyNoise(repIdx)] = pca_MI(MI_opts, R_count_tot, S);
    [MI.nnmf.onlyNoise(repIdx)] = nnmf_MI(MI_opts, R_count_tot, S);

    for shIdx = 1:nShuff
        RSh_count_tot = nan(size(R_count_tot));

        for cellIdx = 1:nNeurons
            RSh_count_tot(cellIdx,:) = R_count_tot(cellIdx,squeeze(shIdxs(shIdx,cellIdx,:)));
        end
        [MISh.linear.onlyNoise(repIdx,shIdx)] = svm_MI(MI_opts, RSh_count_tot', S', cv_type, 'linear_SVM', folds, optimize_params);
        [MISh.RBF.onlyNoise(repIdx,shIdx)] = svm_MI(MI_opts, RSh_count_tot', S', cv_type, 'RBF_SVM', folds, optimize_params);
        [MISh.pca.onlyNoise(repIdx,shIdx)] = pca_MI(MI_opts, RSh_count_tot, S);
        [MISh.nnmf.onlyNoise(repIdx,shIdx)] = nnmf_MI(MI_opts, RSh_count_tot, S);
    end
end

%%
save(['populationInfo_results_',date])
%% Plot figure
load(['populationInfo_results_',date])

maxY = 0.85;

scenario_titles = {'no noise-corr.', 'info-limiting corr.','noise-only info.'};
figure('Position', [360, 347, 785, 271])
for scIdx = 1:numel(scenario_labels)
    scLab = scenario_labels{scIdx};
    count = 0;
    
    subplot(1,3,scIdx)
    hold on
    % Average population information
    meanLin = mean(MI.linear.(scLab));
    errLin = std(MI.linear.(scLab));
    meanLinSh = mean(mean(MISh.linear.(scLab)));
    errLinSh = std(mean(MISh.linear.(scLab),2));
    meanRBF = mean(MI.RBF.(scLab));
    errRBF = std(MI.RBF.(scLab));
    meanRBFSh = mean(mean(MISh.RBF.(scLab)));
    errRBFSh = std(mean(MISh.RBF.(scLab),2));

    meanPCA = mean(MI.pca.(scLab));
    errPCA = std(MI.pca.(scLab));
    meanPCASh = mean(mean(MISh.pca.(scLab)));
    errPCASh = std(mean(MISh.pca.(scLab),2));
    meanNNMF = mean(MI.nnmf.(scLab));
    errNNMF = std(MI.nnmf.(scLab));
    meanNNMFSh = mean(mean(MISh.nnmf.(scLab)));
    errNNMFSh = std(mean(MISh.nnmf.(scLab),2));
    
    % Barplots
    h{1} = bar(scIdx,meanLin,'FaceColor',[0 0 0.7]);
    errorbar(scIdx,meanLin,errLin,'k--')
    h{2} = bar(scIdx+1,meanLinSh,'FaceColor',[0.4 0.4 0.8]);
    errorbar(scIdx+1,meanLinSh,errLinSh,'k--')
    h{3} = bar(scIdx+2,meanRBF,'FaceColor',[0.7 0 0]);
    errorbar(scIdx+2,meanRBF,errRBF,'k--')
    h{4} = bar(scIdx+3,meanRBFSh,'FaceColor',[0.8 0.4 0.4]);
    errorbar(scIdx+3,meanRBFSh,errRBFSh,'k--')
    h{5} = bar(scIdx+4,meanPCA,'FaceColor',[0 0 0.7]);
    errorbar(scIdx+4,meanPCA,errPCA,'k--')
    h{6} = bar(scIdx+5,meanPCASh,'FaceColor',[0.4 0.4 0.8]);
    errorbar(scIdx+5,meanPCASh,errPCASh,'k--')
    h{7} = bar(scIdx+6,meanNNMF,'FaceColor',[0.7 0 0]);
    errorbar(scIdx+6,meanNNMF,errNNMF,'k--')
    h{8} = bar(scIdx+7,meanNNMFSh,'FaceColor',[0.8 0.4 0.4]);
    errorbar(scIdx+7,meanNNMFSh,errNNMFSh,'k--')
    
    ylabel('[bits]')
    if scIdx == 3
        legend([h{1} h{2} h{3} h{4} h{5} h{6} h{7} h{8}], 'Linear', 'Linear shuff.', 'RBF', 'RBF shuff.', 'PCA', 'PCA shuff.', 'NNMF', 'NNMF shuff.')
    end 
    ylim([0,maxY])
    
    title([scenario_titles{scIdx}, ' scenario'])
    
end
sgtitle(['Population information in different correlation scenarios, N = ',num2str(nNeurons),' neurons'])

%% Information breakdown plot
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
    % Compute mean brakdown terms
    meanJoint = mean(mean(MI_breakdown.(scLab).Joint,1),2);
    errJoint = std(squeeze(mean(MI_breakdown.(scLab).Joint,2)));
    meanLin = mean(mean(MI_breakdown.(scLab).ILIN,1),2);
    errLin = std(squeeze(mean(MI_breakdown.(scLab).ILIN,2)));
    meanISS = mean(mean(MI_breakdown.(scLab).ISS,1),2);
    errISS = std(squeeze(mean(MI_breakdown.(scLab).ISS,2)));
    meanICI = mean(mean(MI_breakdown.(scLab).ICI,1),2);
    errICI = std(squeeze(mean(MI_breakdown.(scLab).ICI,2)));
    meanICD = mean(mean(MI_breakdown.(scLab).ICD,1),2);
    errICD = std(squeeze(mean(MI_breakdown.(scLab).ICD,2)));
    
    % Barplots
    % Joint information
    h{1} = bar(scIdx,meanJoint,'FaceColor',[0.4 0.4 0.4]);
    errorbar(scIdx,meanJoint,errJoint,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).Joint,2),[],'tail','right');
    % I_Lin
    h{2} = bar(scIdx+1,meanLin,'FaceColor',[0.8 0.4 0.4]);
    errorbar(scIdx+1,meanLin,errLin,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ILIN,2),[],'tail','right');
    % Sig. sim.
    h{3} = bar(scIdx+2,meanISS,'FaceColor',[0.4 0.8 0.4]);
    errorbar(scIdx+2,meanISS,errISS,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ISS,2),[],'tail','left');
    % Cor-indep
    h{4} = bar(scIdx+3,meanICI,'FaceColor',[0.4 0.4 0.8]);
    errorbar(scIdx+3,meanICI,errICI,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ICI,2));
    % Cor-dep
    h{5} = bar(scIdx+4,meanICD,'FaceColor',[0.8 0.4 0.8]);
    errorbar(scIdx+4,meanICD,errICD,'k--')
    [a,p] = ttest(mean(MI_breakdown.(scLab).ICD,2),[],'tail','right');
    
    ylabel('[bits]')
    if scIdx == 3
        legend([h{1} h{2} h{3} h{4} h{5}], 'Joint', 'Linear', 'Sig. sim.', 'Corr. indep.', 'Corr. dep.')
    end
    ylim([minY,maxY])
    
    title([scenario_titles{scIdx}, ' scenario'])
   
end
sgtitle(['Information breakdown in different correlation scenarios, N = ',num2str(nNeurons),' neurons'])


%% helper functions
function MI = svm_MI(MI_opts, R, S, cv_type, svm_family, folds, optimize_params)
opts.algorithm = svm_family;
opts.cv_type = cv_type;
opts.K = folds;
opts.optimize_params = optimize_params;
opts.verbose = 0; % set to 1 to display asserts and potential Warnings

% Split data into training and test
if strcmp(cv_type,'KFold')
    cvPartition = cvpartition(S,cv_type,opts.K); % CV partitioning used to compute information on the held out trials.
elseif strcmp(cv_type,'LeaveOut')
    cvPartition = cvpartition(S,cv_type,opts.K);
end

% Train and test the SVM
for i = 1:cvPartition.NumTestSets
    test_idxs = find(cvPartition.test(i));
    [tmpPredLabels, labelsTest, dataTest, probs] = buildML(R,S,test_idxs,opts);
    predictedLabels(test_idxs) = tmpPredLabels;
end
info_out = information(predictedLabels, S', MI_opts, {'I'}); % Compute information between real and decoded stimulus, combining the test_sets predicted stimuli to avoid overfitting. We are also bias correcting the information using QE but bias here is very small.
MI = info_out{1};
end

function MI = pca_MI(MIpca_opts, R, S)
reduced_R = pca_pipeline(R, 2, 0);
MI = cell2mat(reduced_R, S, MIpca_opts, {'I'});
end

function MI = nnmf_MI(MIpca_opts, R, S)
reduced_R = nnmf_pipeline(R, 2);
MI = cell2mat(reduced_R, S, MIpca_opts, {'I'});
end