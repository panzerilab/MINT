% -------------------------------------------------------------------------
% Tutorial: Hippocampus_data 
% -------------------------------------------------------------------------
clear; clc
rng(0);
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

for i = 1:1

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
   
    for l = 1:cvPartition.NumTestSets
        test_idxs = find(cvPartition.test(l));
        opts.svm_family = 'linear';
        % Predict Labels Neurons Astrocytes
        [PredLabels_neurons, ~, ~, ~, optimParams] = svm_pipeline(double(R_neuron)',S',test_idxs, opts);
        [PredLabels_astro, ~, ~, ~, optimParams] = svm_pipeline(double(R_astro)',S',test_idxs, opts);
        % Save Labels
        predLabelsLin.(field_name).unshuffled.neurons(l) = {PredLabels_neurons'};
        predLabelsLin.(field_name).unshuffled.astros(l)  = {PredLabels_astro'};
        predLabelsLin.(field_name).unshuffled.s(l)       = {S(test_idxs)};
        predictedLabels = [PredLabels_neurons, PredLabels_astro];

        % Calculate Information
        info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
        info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
        info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});

        info_out_all_n(end+1) = info_out_neuron{1};
        info_out_all_a(end+1) = info_out_astro{1};
        info_out_all_na(end+1) = info_out_na{1};
    end
    info_out_mean_n = mean(info_out_all_n);
    disp(info_out_mean_n)
    info_out_mean_a = mean(info_out_all_a);
    info_out_mean_na = mean(info_out_all_na);

    MI_n.linear.(field_name).mean = info_out_mean_n;
    MI_n.linear.(field_name).all = info_out_all_n;

    MI_a.linear.(field_name).mean = info_out_mean_a;
    MI_a.linear.(field_name).all = info_out_all_a;

    MI_na.linear.(field_name).mean = info_out_mean_na;
    MI_na.linear.(field_name).all = info_out_all_na;
    
    %RBF SVM
    info_out_all_n = [];
    info_out_all_a = [];
    info_out_all_na = [];

    for r = 1:cvPartition.NumTestSets
        test_idxs = find(cvPartition.test(r));
        opts.svm_family = 'RBF';
        [PredLabels_neurons, ~, ~, ~, optimParams] = svm_pipeline(double(R_neuron)',S',test_idxs, opts);
        [PredLabels_astro, ~, ~, ~, optimParams] = svm_pipeline(double(R_astro)',S',test_idxs, opts);

        predLabelsRBF.(field_name).unshuffled.neurons(r) = {PredLabels_neurons'};
        predLabelsRBF.(field_name).unshuffled.astros(r)  = {PredLabels_astro'};
        predLabelsRBF.(field_name).unshuffled.s(r)       = {S(test_idxs)};
        predictedLabels = [PredLabels_neurons, PredLabels_astro];

        info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
        info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
        info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});

        info_out_all_n(end+1) = info_out_neuron{1};
        info_out_all_a(end+1) = info_out_astro{1};
        info_out_all_na(end+1) = info_out_na{1};
    end
    info_out_mean_n = mean(info_out_all_n);
    info_out_mean_a = mean(info_out_all_a);
    info_out_mean_na = mean(info_out_all_na);

    MI_n.RBF.(field_name).mean = info_out_mean_n;
    MI_n.RBF.(field_name).all = info_out_all_n;

    MI_a.RBF.(field_name).mean = info_out_mean_a;
    MI_a.RBF.(field_name).all = info_out_all_a;

    MI_na.RBF.(field_name).mean = info_out_mean_na;
    MI_na.RBF.(field_name).all = info_out_all_na;

    % Perform shufflung and compute MI for shuffled data
    for shIdx = 1:nShuff
        RSh_neuron = shuffle(S', R_neuron', 1, [1,0]);
        RSh_astro = shuffle(S', R_astro', 1, [1,0]);

        % Linear SVM
        info_out_all_n = [];
        info_out_all_a = [];
        info_out_all_na = [];
        for l = 1:cvPartition.NumTestSets
            test_idxs = find(cvPartition.test(l));
            opts.svm_family = 'linear';
            [PredLabels_neurons, ~ , ~ , ~,optimParams ] = svm_pipeline(double(RSh_neuron),S',test_idxs, opts);
            [PredLabels_astro, ~ , ~ , ~, ~ ] = svm_pipeline(double(RSh_astro),S',test_idxs, opts);

            predLabelsLin.(field_name).shuffled.neurons(l) = {PredLabels_neurons'};
            predLabelsLin.(field_name).shuffled.astros(l)  = {PredLabels_astro'};
            predLabelsLin.(field_name).shuffled.s(l)       = {S(test_idxs)};
            predictedLabels = [PredLabels_neurons, PredLabels_astro];

            info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
            info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
            info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});
    
            info_out_all_n(end+1) = info_out_neuron{1};
            info_out_all_a(end+1) = info_out_astro{1};
            info_out_all_na(end+1) = info_out_na{1};
        end
        info_out_mean_n = mean(info_out_all_n);
        info_out_mean_a = mean(info_out_all_a);
        info_out_mean_na = mean(info_out_all_na);

        MISh_n.linear.(field_name)(shIdx).mean = info_out_mean_n;
        MISh_n.linear.(field_name)(shIdx).all = info_out_all_n;
        
        MISh_a.linear.(field_name)(shIdx).mean = info_out_mean_a;
        MISh_a.linear.(field_name)(shIdx).all = info_out_all_a;

        MISh_na.linear.(field_name)(shIdx).mean = info_out_mean_na;
        MISh_na.linear.(field_name)(shIdx).all = info_out_all_na;

        % RBF SVM
        info_out_all_n = [];
        info_out_all_a = [];
        info_out_all_na = [];

        for r = 1:cvPartition.NumTestSets
            test_idxs = find(cvPartition.test(r));
            opts.svm_family = 'RBF';
            [PredLabels_neurons, ~ , ~ , ~,optimParams ] = svm_pipeline(double(RSh_neuron),S',test_idxs, opts);
            [PredLabels_astro, ~ , ~ , ~, ~ ] = svm_pipeline(double(RSh_astro),S',test_idxs, opts);
            
            predLabelsRBF.(field_name).shuffled.neurons(r) = {PredLabels_neurons'};
            predLabelsRBF.(field_name).shuffled.astros(r)  = {PredLabels_astro'};
            predLabelsRBF.(field_name).shuffled.s(r)       = {S(test_idxs)};
            predictedLabels = [PredLabels_neurons, PredLabels_astro];

            info_out_neuron = information(PredLabels_neurons', S(test_idxs), MI_opts, {'I'});
            info_out_astro = information(PredLabels_astro', S(test_idxs), MI_opts, {'I'});
            info_out_na = information(predictedLabels', S(test_idxs), MI_opts, {'I'});

            info_out_all_n(end+1) = info_out_neuron{1};
            info_out_all_a(end+1) = info_out_astro{1};
            info_out_all_na(end+1) = info_out_na{1};
        end
        info_out_mean_n = mean(info_out_all_n);
        info_out_mean_a = mean(info_out_all_a);
        info_out_mean_na = mean(info_out_all_na);

        MISh_n.RBF.(field_name)(shIdx).mean = info_out_mean_n;
        MISh_n.RBF.(field_name)(shIdx).all = info_out_all_n;
        
        MISh_a.RBF.(field_name)(shIdx).mean = info_out_mean_a;
        MISh_a.RBF.(field_name)(shIdx).all = info_out_all_a;

        MISh_na.RBF.(field_name)(shIdx).mean = info_out_mean_na;
        MISh_na.RBF.(field_name)(shIdx).all = info_out_all_na;
    end
end


%%


% % Calculate Mean Information Breakdown
% meanJoint_all = [];
% meanLin_all = [];
% meanISS_all = [];
% meanICI_all = [];
% meanICD_all = [];
% 
% titles = field_names';
% 
% for subIdx = 1:numel(field_names)
%     subLab = field_names{subIdx,1};
%     meanJoint_all = [meanJoint_all, mean(mean(MI_breakdown.(subLab).Joint,1),2)];
%     meanLin_all = [meanLin_all, mean(mean(MI_breakdown.(subLab).ILIN,1),2)];
%     meanISS_all = [meanISS_all, mean(mean(MI_breakdown.(subLab).ISS,1),2)];
%     meanICI_all = [meanICI_all, mean(mean(MI_breakdown.(subLab).ICI,1),2)];
%     meanICD_all = [meanICD_all, mean(mean(MI_breakdown.(subLab).ICD,1),2)];
% end
% %%
% [meanJoint, errJoint] = compute_stats(meanJoint_all);
% [meanLin, errLin] = compute_stats(meanLin_all);
% [meanISS, errISS] = compute_stats(meanISS_all);
% [meanICI, errICI] = compute_stats(meanICI_all);
% [meanICD, errICD] = compute_stats(meanICD_all);

%%
% Calculate Mean Population Information 
% MILinear_all = [];
% MILinearSh_all = [];
% MIRBF_all = [];
% MIRBFSh_all = [];
% 
% for subIdx = 1:numel(field_names)
%     subLab = field_names{subIdx,1};
%     MILinear_all = [MILinear_all, mean(MI.linear.(subLab))];
%     MILinearSh_all = [MILinearSh_all, mean(MISh.linear.(subLab))];
%     MIRBF_all = [MIRBF_all, mean(MI.RBF.(subLab))];
%     MIRBFSh_all = [MIRBFSh_all, mean(MISh.RBF.(subLab))];
% end
% %%
% [meanMILin, errMILin] = compute_stats(MILinear_all);
% [meanMILinSh, errMILinSh] = compute_stats(MILinearSh_all);
% [meanMIRBF, errMIRBF] = compute_stats(MIRBF_all);
% [meanMIRBFSh, errMIRBFSh] = compute_stats(MIRBFSh_all);
% %%
% % save('meanLin_astro.mat', 'meanLin_all');
% % save('meanJoint_astro.mat', 'meanJoint_all');
% % save('meanISS_astro.mat', 'meanISS_all');
% % save('meanICI_astro.mat', 'meanICI_all');
% % save('meanICD_astro.mat', 'meanICD_all');
% save('MILinear_astro.mat', 'MILinear_all');
% save('MILinearSh_astro.mat', 'MILinearSh_all');
% save('MIRBF_astro.mat', 'MIRBF_all');
% save('MIRBFSh_astro.mat', 'MIRBFSh_all');
% 
% %% Plot 
% fSize = 15;
% 
% % % Plot Information Breakdown
% % maxY = 0.03;
% % minY = -0.005;
% % 
% % [~, p_ISS] = ttest(meanISS_all, 0);
% % [~, p_ICI] = ttest(meanICI_all, 0);
% % [~, p_ICD] = ttest(meanICD_all, 0);
% % 
% % figure('Position', [360, 347, 785, 271])
% % hold on 
% % % Barplots
% % % Joint information
% % h{1} = bar(1, meanJoint,'FaceColor',[0.4 0.4 0.4]);
% % errorbar(1, meanJoint,errJoint,'k--')
% % % I_Lin
% % h{2} = bar(2, meanLin,'FaceColor',[0.8 0.4 0.4]);
% % errorbar(2, meanLin,errLin,'k--')
% % % Sig. sim.
% % h{3} = bar(3, meanISS,'FaceColor',[0.4 0.8 0.4]);
% % errorbar(3, meanISS,errISS,'k--')
% % % Cor-indep
% % h{4} = bar(4, meanICI,'FaceColor',[0.4 0.4 0.8]);
% % errorbar(4, meanICI,errICI,'k--');
% % % Cor-dep
% % h{5} = bar(5, meanICD,'FaceColor',[0.8 0.4 0.8]);
% % errorbar(5, meanICD,errICD,'k--');
% % 
% % xticks([]);
% % set(gca, 'FontSize', 12);
% % ylabel('Information [bits]')
% % legend_obj = legend([h{1} h{2} h{3} h{4} h{5}], 'Joint', 'Linear', 'Sig. sim.', 'Corr. indep.', 'Corr. dep.');
% % legend_obj.FontSize = 14;
% % legend('AutoUpdate','off');
% % ylim([minY, maxY])
% % set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 60]);
% % ax = gca;
% % ax.FontSize = 12;
% % title('Information Breakdown Astrocytes', 'FontSize', fSize+5);
% % 
% % add_pvalue([p_ISS, p_ICI, p_ICD], maxY, [2.75,3.75,4.75], [3.25,4.25,5.25], 0.0005, [0.025, 0.025, 0.02])
% % 
% % saveas(gcf, 'HippocampusData_Astrocytes_InformationBreakdown.png');
% 
% % Plot Population Information
% maxY = 0.7;
% minY = 0;
% 
% [~, p_LinRBF] = ttest(MILinear_all, MIRBF_all);
% [~, p_RBFsh] = ttest(MIRBF_all, MIRBFSh_all);
% [~, p_Linsh] = ttest(MILinear_all, MILinearSh_all);
% 
% 
% figure('Position', [360, 347, 785, 271])
% hold on
% % Barplots
% h{1} = bar(1,meanMILin,'FaceColor',[0 0 0.7]);
% errorbar(1,meanMILin,errMILin,'k--', 'LineWidth', 1)
% h{2} = bar(2,meanMILinSh,'FaceColor',[0.4 0.4 0.8]);
% errorbar(2,meanMILinSh,errMILinSh,'k--', 'LineWidth', 1)
% h{3} = bar(3,meanMIRBF,'FaceColor',[0.7 0 0]);
% errorbar(3,meanMIRBF,errMIRBF,'k--', 'LineWidth', 1)
% h{4} = bar(4,meanMIRBFSh,'FaceColor',[0.8 0.4 0.4]);
% errorbar(4,meanMIRBFSh,errMIRBFSh,'k--', 'LineWidth', 1)
% 
% 
% xticks([]);
% set(gca, 'FontSize', 12);
% ylabel('Information [bits]')
% legend_obj = legend([h{1} h{2} h{3} h{4}], 'Linear', 'Linear shuff.', 'RBF', 'RBF shuff.');
% legend_obj.FontSize = 14;
% legend('AutoUpdate','off');
% 
% ylim([minY, maxY])
% set(gcf, 'Units', 'centimeters', 'Position', [0, 0, 40, 60]);
% ax = gca;
% ax.FontSize = 12;
% 
% title_text = title('Population Information Astrocytes', 'FontSize', fSize + 5);
% title_position = get(title_text, 'Position');
% title_position(2) = maxY+0.01;  
% set(title_text, 'Position', title_position);
% 
% add_pvalue([p_LinRBF, p_Linsh, p_RBFsh], maxY, [1,1,3], [3,2,4], 0.005, [0.1, 0.2, 0.15])
% saveas(gcf, 'HippocampusData_PopulationInformation_Astrocytes.png', 'd');
% 
% %% Helper Function
% 
% function [mean_val, sem_val] = compute_stats(data)
%     mean_val = mean(data);
%     std_val = std(data);
%     n = length(data);
%     sem_val = std_val / sqrt(n);
% end
% 
% function add_pvalue(p_values, y, x1, x2, text_heights, line_heights)
%     for i = 1:numel(p_values)
%         stars = get_stars(p_values(i));
%         line([x1(i), x2(i)], [y - line_heights(i), y - line_heights(i)], 'Color', 'black', 'LineWidth', 1);
%         text((x1(i) + x2(i)) / 2, y - line_heights(i) + text_heights, stars, 'HorizontalAlignment', 'center', 'Color', 'black');
%     end
% end
% 
% function stars = get_stars(p_value)
%     if p_value < 0.001
%         stars = '***';
%     elseif p_value < 0.01
%         stars = '**';
%     elseif p_value < 0.05
%         stars = '*';
%     else
%         stars = 'n.s.';
%     end
% end
% 
