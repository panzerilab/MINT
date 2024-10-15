
% Tutorial: Hippocampus_data 
% -------------------------------------------------------------------------
clear; clc
% rng(0);
disp('Astrocyte Data')
load('Hippocampus_data.mat');

field_names = fieldnames(Hippocampus_data);

% Population MI options (used for SVM)
MI_opts.verbose = false;
MI_opts.method = "dr";
MI_opts.bias = 'qe';
MI_opts.btsp = 0; 
MI_opts.bin_methodX = 'none'; 
MI_opts.bin_methodY = 'none';

% Info breakdown opts
% info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
nShuff = 2;
predLabelsLin = {};
predLabelsRBF = {};



for i=7:7 %1:length(field_names)
    % -------------------------------------------------------------------------
    % Step 1: Load Subject data
    % -------------------------------------------------------------------------
    
    % Extract data
    field_name = field_names{i, 1};
    disp(['Subject ',field_name])
    S = Hippocampus_data.(field_name).Position;
    R_neuron = Hippocampus_data.(field_name).R_neuron;
    R_neuron_binned = Hippocampus_data.(field_name).R_neuron_binned;
    R_astro = Hippocampus_data.(field_name).R_astro;
    R_astro_binned = Hippocampus_data.(field_name).R_astro;
    
    [nAstro, ~] = size(R_astro);
    [nNeurons, nTrials] = size(R_neuron);
    nPairs = nAstro*(nAstro-1)/2;
    
    % %Initialize Structures
    % for bdwIdx = 1:numel(info_bdw_terms)
    %     bdwLab = info_bdw_terms{bdwIdx};
    %     MI_breakdown.(field_name).(bdwLab) = nan(1, nPairs);
    % end
    % noiseCorr.(field_name) = nan(1, 2, nPairs);
    
    % -------------------------------------------------------------------------
    % Step 2: Compute Noise Correlation and Information Breakdown
    % -------------------------------------------------------------------------
    
    % pairIdx = 0;
    % for cell1 = 1:nAstro
    %     for cell2 = cell1+1:nAstro
    %         pairIdx = pairIdx+1;
    %         jointResp = [R_binned(cell1,:);R_binned(cell2,:)];
    %         if ~isempty(cell2)
    %             infoBdw = information(jointResp,S,MI_opts,{'I','ILIN','ISS','ICI','ICD'});
    %         end
    % 
    %         % Store information breakdown results
    %         for bdwIdx = 1:numel(info_bdw_terms)
    %             bdwLab = info_bdw_terms{bdwIdx};
    %             MI_breakdown.(field_name).(bdwLab)(pairIdx) = infoBdw{bdwIdx};
    %         end
    %     end
    % end
    
    % -------------------------------------------------------------------------
    % Step 3: Compute Population Information with SVM Decoder
    % -------------------------------------------------------------------------
    
    
    
    opts.optimize_params = true;
    opts.parallel_optim = true;
    % if opts.parallel_optim == true
    %     if isempty(gcp('nocreate'))
    %        parpool('local', 8);
    %     end
    % end 
    opts.cv_type = "KFold"; 
    opts.K = 5; 
    opts.libsvm = false; 
    
    % Split data into training and test
    cvPartition = cvpartition(S, 'KFold', 10); 
    
    %Linear SVM
    info_out_all_n = [];
    info_out_all_a = [];
    info_out_all_na = [];
    info_out_all_nap = [];
    
    
    % for l = 1:cvPartition.NumTestSets
    %     test_idxs = find(cvPartition.test(l));
    %     opts.svm_family = 'linear';
    %     % Predict Labels Neurons Astrocytes
    %     [PredLabels_pooled, ~, ~, ~, optimParams] = svm_pipeline(double([R_neuron; R_astro])',S',test_idxs, opts);
    %     [PredLabels_neurons, ~, ~, ~, optimParams] = svm_pipeline(double(R_neuron)',S',test_idxs, opts);
    %     [PredLabels_astro, ~, ~, ~, optimParams] = svm_pipeline(double(R_astro)',S',test_idxs, opts);
    %     % Save Labels
    %     predLabelsLin.(field_name).unshuffled.neurons(l) = {PredLabels_neurons'};
    %     predLabelsLin.(field_name).unshuffled.astros(l)  = {PredLabels_astro'};
    %     predLabelsLin.(field_name).unshuffled.s(l)       = {S(test_idxs)};
    %     predictedLabels = [PredLabels_neurons, PredLabels_astro];
    % 
    %     % Calculate Information
    %     info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
    %     info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
    %     info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});
    %     info_out_nap = information(PredLabels_pooled', S(test_idxs), MI_opts, {'I'});
    % 
    % 
    %     info_out_all_n(end+1) = info_out_neuron{1};
    %     info_out_all_a(end+1) = info_out_astro{1};
    %     info_out_all_na(end+1) = info_out_na{1};
    %     info_out_all_nap(end+1) = info_out_nap{1};
    % end
    % info_out_mean_n = mean(info_out_all_n);
    % info_out_mean_a = mean(info_out_all_a);
    % info_out_mean_na = mean(info_out_all_na);
    % info_out_mean_nap = mean(info_out_all_nap);
    % 
    % 
    % MI_n.linear.(field_name).mean = info_out_mean_n;
    % MI_n.linear.(field_name).all = info_out_all_n;
    % 
    % MI_a.linear.(field_name).mean = info_out_mean_a;
    % MI_a.linear.(field_name).all = info_out_all_a;
    % 
    % MI_na.linear.(field_name).mean = info_out_mean_na;
    % MI_na.linear.(field_name).all = info_out_all_na;
    % 
    % MI_nap.linear.(field_name).mean = info_out_mean_nap;
    % MI_nap.linear.(field_name).all = info_out_all_nap;
    % 
    % %RBF SVM
    % info_out_all_n = [];
    % info_out_all_a = [];
    % info_out_all_na = [];
    % info_out_all_nap = [];
    % 
    % for r = 1:cvPartition.NumTestSets
    %     test_idxs = find(cvPartition.test(r));
    %     opts.svm_family = 'RBF';
    %     [PredLabels_pooled, ~, ~, ~, optimParams] = svm_pipeline(double([R_neuron; R_astro])',S',test_idxs, opts);
    %     [PredLabels_neurons, ~, ~, ~, optimParams] = svm_pipeline(double(R_neuron)',S',test_idxs, opts);
    %     [PredLabels_astro, ~, ~, ~, optimParams] = svm_pipeline(double(R_astro)',S',test_idxs, opts);
    % 
    %     predLabelsRBF.(field_name).unshuffled.neurons(r) = {PredLabels_neurons'};
    %     predLabelsRBF.(field_name).unshuffled.astros(r)  = {PredLabels_astro'};
    %     predLabelsRBF.(field_name).unshuffled.s(r)       = {S(test_idxs)};
    %     predictedLabels = [PredLabels_neurons, PredLabels_astro];
    % 
    %     info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
    %     info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
    %     info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});
    %     info_out_nap = information(PredLabels_pooled', S(test_idxs), MI_opts, {'I'});
    % 
    %     info_out_all_n(end+1) = info_out_neuron{1};
    %     info_out_all_a(end+1) = info_out_astro{1};
    %     info_out_all_na(end+1) = info_out_na{1};
    %     info_out_all_nap(end+1) = info_out_nap{1};
    % end
    % info_out_mean_n = mean(info_out_all_n);
    % info_out_mean_a = mean(info_out_all_a);
    % info_out_mean_na = mean(info_out_all_na);
    % info_out_mean_nap = mean(info_out_all_nap);
    % 
    % MI_n.RBF.(field_name).mean = info_out_mean_n;
    % MI_n.RBF.(field_name).all = info_out_all_n;
    % 
    % MI_a.RBF.(field_name).mean = info_out_mean_a;
    % MI_a.RBF.(field_name).all = info_out_all_a;
    % 
    % MI_na.RBF.(field_name).mean = info_out_mean_na;
    % MI_na.RBF.(field_name).all = info_out_all_na;
    % 
    % MI_nap.RBF.(field_name).mean = info_out_mean_nap;
    % MI_nap.RBF.(field_name).all = info_out_all_nap;
    % 
    
    % Perform shufflung and compute MI for shuffled data
    for shIdx = 1:nShuff
        RSh_neuron = shuffle(S', R_neuron', 1, [1,0]);
        RSh_astro = shuffle(S', R_astro', 1, [1,0]);
    
        % Linear SVM
        info_out_all_n = [];
        info_out_all_a = [];
        info_out_all_na = [];
        info_out_all_nap = [];
        for l = 1:cvPartition.NumTestSets
            test_idxs = find(cvPartition.test(l));
            opts.svm_family = 'linear';
            [PredLabels_pooled, ~, ~, ~, optimParams] = svm_pipeline(double([RSh_neuron, RSh_astro]),S',test_idxs, opts);
            [PredLabels_neurons, ~ , ~ , ~,optimParams ] = svm_pipeline(double(RSh_neuron),S',test_idxs, opts);
            [PredLabels_astro, ~ , ~ , ~, ~ ] = svm_pipeline(double(RSh_astro),S',test_idxs, opts);
    
            predLabelsLin.(field_name).shuffled.neurons(l) = {PredLabels_neurons'};
            predLabelsLin.(field_name).shuffled.astros(l)  = {PredLabels_astro'};
            predLabelsLin.(field_name).shuffled.s(l)       = {S(test_idxs)};
            predictedLabels = [PredLabels_neurons, PredLabels_astro];
    
            info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
            info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
            info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});
            info_out_nap = information(PredLabels_pooled', S(test_idxs), MI_opts, {'I'});
    
            info_out_all_n(end+1) = info_out_neuron{1};
            info_out_all_a(end+1) = info_out_astro{1};
            info_out_all_na(end+1) = info_out_na{1};
            info_out_all_nap(end+1) = info_out_nap{1};
        end
        info_out_mean_n = mean(info_out_all_n);
        info_out_mean_a = mean(info_out_all_a);
        info_out_mean_na = mean(info_out_all_na);
        info_out_mean_nap = mean(info_out_all_nap);
    
    
        MISh_n.linear.(field_name)(shIdx).mean = info_out_mean_n;
        MISh_n.linear.(field_name)(shIdx).all = info_out_all_n;
        
        MISh_a.linear.(field_name)(shIdx).mean = info_out_mean_a;
        MISh_a.linear.(field_name)(shIdx).all = info_out_all_a;
    
        MISh_na.linear.(field_name)(shIdx).mean = info_out_mean_na;
        MISh_na.linear.(field_name)(shIdx).all = info_out_all_na;
    
        MISh_nap.linear.(field_name)(shIdx).mean = info_out_mean_nap;
        MISh_nap.linear.(field_name)(shIdx).all = info_out_all_nap;
    
        % RBF SVM
        info_out_all_n = [];
        info_out_all_a = [];
        info_out_all_na = [];
        info_out_all_nap = [];
        for r = 1:cvPartition.NumTestSets
            test_idxs = find(cvPartition.test(r));
            opts.svm_family = 'RBF';
            [PredLabels_pooled, ~, ~, ~, optimParams] = svm_pipeline(double([RSh_neuron, RSh_astro]),S',test_idxs, opts);
            [PredLabels_neurons, ~ , ~ , ~,optimParams ] = svm_pipeline(double(RSh_neuron),S',test_idxs, opts);
            [PredLabels_astro, ~ , ~ , ~, ~ ] = svm_pipeline(double(RSh_astro),S',test_idxs, opts);
            
            predLabelsRBF.(field_name).shuffled.neurons(r) = {PredLabels_neurons'};
            predLabelsRBF.(field_name).shuffled.astros(r)  = {PredLabels_astro'};
            predLabelsRBF.(field_name).shuffled.s(r)       = {S(test_idxs)};
            predictedLabels = [PredLabels_neurons, PredLabels_astro];
    
            info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
            info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
            info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});
            info_out_nap = information(PredLabels_pooled', S(test_idxs), MI_opts, {'I'});
            
            info_out_all_n(end+1) = info_out_neuron{1};
            info_out_all_a(end+1) = info_out_astro{1};
            info_out_all_na(end+1) = info_out_na{1};
            info_out_all_nap(end+1) = info_out_nap{1};
        end
        info_out_mean_n = mean(info_out_all_n);
        info_out_mean_a = mean(info_out_all_a);
        info_out_mean_na = mean(info_out_all_na);
        info_out_mean_nap = mean(info_out_all_nap);
    
        MISh_n.RBF.(field_name)(shIdx).mean = info_out_mean_n;
        MISh_n.RBF.(field_name)(shIdx).all = info_out_all_n;
        
        MISh_a.RBF.(field_name)(shIdx).mean = info_out_mean_a;
        MISh_a.RBF.(field_name)(shIdx).all = info_out_all_a;
    
        MISh_na.RBF.(field_name)(shIdx).mean = info_out_mean_na;
        MISh_na.RBF.(field_name)(shIdx).all = info_out_all_na;
    
        MISh_nap.RBF.(field_name)(shIdx).mean = info_out_mean_nap;
        MISh_nap.RBF.(field_name)(shIdx).all = info_out_all_nap;
    end
end
%%
%save('MI_n.mat', 'MI_n')
%save('MI_a.mat', 'MI_a')
%save('MI_na.mat', 'MI_na')
%save('MI_nap.mat', 'MI_nap')

%save('MISh_n.mat', 'MISh_n')
%save('MISh_a.mat', 'MISh_a')
%save('MISh_na.mat', 'MISh_na')
%save('MISh_nap.mat', 'MISh_nap')

save('predLabelsLin.mat', 'predLabelsLin')
save('predLabelsRBF.mat', 'predLabelsRBF')

