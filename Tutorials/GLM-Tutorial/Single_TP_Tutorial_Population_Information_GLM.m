clc, clear,warning off,
rng(123);
total_time_start = tic;

% Load Data
load('AC_data.mat');
load('PPC_data.mat');

datasetsAC = length(AC_data);
datasetsPPC = length(PPC_data);
[~,~,nTP_choice] = size(AC_data(1).neural_data_choice);
[~,~,nTP_stimulus] = size(AC_data(1).neural_data_stimulus);

% Initialize structures
info_out_choice_AC = [];
info_out_choice_PPC = [];
info_out_stimulus_AC = [];
info_out_stimulus_PPC = [];
info_out_choice_AC_all = zeros(datasetsAC, nTP_choice);
info_out_choice_PPC_all = zeros(datasetsPPC,nTP_choice);
info_out_stimulus_AC_all = zeros(datasetsAC, nTP_stimulus);
info_out_stimulus_PPC_all = zeros(datasetsPPC,nTP_stimulus);

% Input Parameter for the GLM Model 
opts.distribution = 'binomial';
opts.alpha = 0.95;
opts.lambda = [0.001, 0.005, 0.01, 0.05, 0.1, 0.5, 1];
opts.CV = false;
kFolds = 5;
opts.regularization = 'elasticNet';
opts.glmnet=true;

% MI options
MI_opts = struct('verbose', false, 'method', 'dr', 'bias', 'pt');

% -------------------------------------------------------------------------
% 1. Stimulus and Choice Information in AC
% -------------------------------------------------------------------------

disp('Calculating Stimulus and Choice Information in AC');

for dataset = 1:length(AC_data)
    name = AC_data(dataset).data;
    disp([num2str(dataset),'.Dataset: ', name]);
    
    % load dataset data
    neural_data_choice = AC_data(dataset).neural_data_choice;
    neural_data_stimulus = AC_data(dataset).neural_data_stimulus;
    choice = AC_data(dataset).choice;
    stimulus = AC_data(dataset).stimulus;

    % Initialize structures and variables
    [numTrials, numCells, numTimePoints_choice] = size(neural_data_choice);
    [~,~,numTimePoints_stimulus] = size(neural_data_stimulus);
     predicted_choice_AC = zeros(numTimePoints_choice, numTrials);
    predicted_stimulus_AC = zeros(numTimePoints_stimulus, numTrials);

    % Predict choice for each frame 
    for tPc = 1:numTimePoints_choice
        for i = 1:length(fields(AC_data(dataset).trainIndices))
            fold_idx = ['Fold_', num2str(i)];
            test_idxs = AC_data(dataset).testIndices.(fold_idx);
            [tmpPredLabels] = glm_pipeline(neural_data_choice(:, :, tPc),choice, test_idxs, opts);
            predicted_choice_AC(tPc, test_idxs)= tmpPredLabels;
        end
    end  
    predicted_choice_AC = predicted_choice_AC';
   
    % Predict stimulus for each frame 
    for tPs = 1:numTimePoints_stimulus
        for i = 1:length(fields(AC_data(dataset).trainIndices))
            fold_idx = ['Fold_', num2str(i)];
            test_idxs = AC_data(dataset).testIndices.(fold_idx);
            [tmpPredLabels] = glm_pipeline(neural_data_stimulus(:, :, tPs), stimulus, test_idxs, opts);
            predicted_stimulus_AC(tPs, test_idxs)= tmpPredLabels;
        end
    end
    predicted_stimulus_AC = predicted_stimulus_AC';

    % Compute information between real and decoded stimulus
    for tPc = 1:numTimePoints_choice
        predicted_choice_tmp = double(mean(predicted_choice_AC(:, tPc), 2) >= 0.5);
        info_out = information(predicted_choice_tmp', choice', MI_opts, {'I'});
        info_out_choice_AC(tPc) =  info_out{1, 1};
    end
    info_out_choice_AC_all(dataset,:) = info_out_choice_AC;

    % Compute information between real and decoded stimulus
    for tPs = 1:numTimePoints_stimulus
        predicted_stimulus_tmp = double(mean(predicted_stimulus_AC(:, tPs), 2) >= 0.5);
        info_out = information(predicted_stimulus_tmp', stimulus', MI_opts, {'I'});
        info_out_stimulus_AC(tPs) =  info_out{1, 1};
    end
    info_out_stimulus_AC_all(dataset,:) = info_out_stimulus_AC;
end

% -------------------------------------------------------------------------
% 2. Stimulus and Choice Information in PPC
% -------------------------------------------------------------------------

disp('Calculating Stimulus and Choice Information in PPC');

for dataset = 1:length(PPC_data)
    name = PPC_data(dataset).data;
    disp([num2str(dataset),'.Dataset: ', name]);

    % load dataset data
    neural_data_choice = PPC_data(dataset).neural_data_choice;
    neural_data_stimulus = PPC_data(dataset).neural_data_stimulus;
    choice = PPC_data(dataset).choice;
    stimulus = PPC_data(dataset).stimulus;
    
    % Initialize structures and variables
    [numTrials, numCells, numTimePoints_choice] = size(neural_data_choice);
    [~,~,numTimePoints_stimulus] = size(neural_data_stimulus);
    predicted_choice_PPC = zeros(numTimePoints_choice, numTrials);
    predicted_stimulus_PPC = zeros(numTimePoints_stimulus, numTrials);
    
    % Predict choice for each frame
    for tPc = 1:numTimePoints_choice
        for i = 1:length(fields(PPC_data(dataset).trainIndices))
            fold_idx = ['Fold_', num2str(i)];
            test_idxs = PPC_data(dataset).testIndices.(fold_idx);
            [tmpPredLabels] = glm_pipeline(neural_data_choice(:, :, tPc), choice, test_idxs, opts);
            predicted_choice_PPC(tPc, test_idxs)= tmpPredLabels;
        end
    end
    predicted_choice_PPC = predicted_choice_PPC';
    
    % Predict stimulus for each frame
    for tPs = 1:numTimePoints_stimulus
        for i = 1:length(fields(PPC_data(dataset).trainIndices))
            fold_idx = ['Fold_', num2str(i)];
            test_idxs = PPC_data(dataset).testIndices.(fold_idx);
            [tmpPredLabels] = glm_pipeline(neural_data_stimulus(:, :, tPs), stimulus, test_idxs, opts);
            predicted_stimulus_PPC(tPs, test_idxs)= tmpPredLabels;
        end
    end
    predicted_stimulus_PPC = predicted_stimulus_PPC';
    

    % Compute information between real and decoded choice
    for t = 1:numTimePoints_choice
        predicted_choice_tmp = double(mean(predicted_choice_PPC(:, t), 2) >= 0.5);
        info_out = information(predicted_choice_tmp', choice', MI_opts, {'I'});
        info_out_choice_PPC(t) =  info_out{1, 1};
    end
    info_out_choice_PPC_all(dataset,:) = info_out_choice_PPC;

    % Compute information between real and decoded stimulus
    for t = 1:numTimePoints_stimulus
        predicted_stimulus_tmp = double(mean(predicted_stimulus_PPC(:, t), 2) >= 0.5);
        info_out = information(predicted_stimulus_tmp', stimulus', MI_opts, {'I'});
        info_out_stimulus_PPC(t) =  info_out{1, 1};
    end
    info_out_stimulus_PPC_all(dataset,:) = info_out_stimulus_PPC;
end

% Calculate Mean Information of all Datasets and SEM 
[mean_AC_choice_information, sem_AC_choice] = compute_stats(info_out_choice_AC_all);
[mean_AC_stimulus_information, sem_AC_stimulus] = compute_stats(info_out_stimulus_AC_all);

[mean_PPC_choice_information, sem_PPC_choice] = compute_stats(info_out_choice_PPC_all);
[mean_PPC_stimulus_information, sem_PPC_stimulus] = compute_stats(info_out_stimulus_PPC_all);


total_time_end = toc(total_time_start);
fprintf('Total Calculation Time =: %.4f Sekunden\n', total_time_end);

%% Plot
samplingRate = 15.6;
x_values_choice = ((1:length(mean_PPC_choice_information))/15.6);
x_values_stimulus = ((1:length(mean_PPC_stimulus_information))/15.6);

% Plot choice information above time 
figure;
subplot(2, 1, 1);
hold on;
p1 = plot(x_values_choice, mean_PPC_choice_information, 'LineWidth', 2, 'DisplayName', 'PPC');
p2 = plot(x_values_choice, mean_AC_choice_information, 'LineWidth', 2, 'DisplayName', 'AC');
fill([x_values_choice, fliplr(x_values_choice)], [mean_PPC_choice_information - sem_PPC_choice, ...
    fliplr(mean_PPC_choice_information + sem_PPC_choice)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([x_values_choice, fliplr(x_values_choice)], [mean_AC_choice_information - sem_AC_choice, ...
    fliplr(mean_AC_choice_information + sem_AC_choice)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xline(47/samplingRate, '--', 'LineWidth', 1.5, 'Color', 'k');
hold off;
xlabel('Time (s)');
ylabel('Information (bits)');
title('Cumulative Choice Category Information');
legend([p1, p2], 'Location', 'Best');
xlim([0, 6.2]);

% Plot stimulus information above time
subplot(2, 1, 2);
hold on;
p3 = plot(x_values_stimulus, mean_PPC_stimulus_information, 'LineWidth', 2, 'DisplayName', 'PPC');
p4 = plot(x_values_stimulus, mean_AC_stimulus_information, 'LineWidth', 2, 'DisplayName', 'AC');
fill([x_values_stimulus, fliplr(x_values_stimulus)], [mean_PPC_stimulus_information - sem_PPC_stimulus, ...
    fliplr(mean_PPC_stimulus_information + sem_PPC_stimulus)], 'b', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
fill([x_values_stimulus, fliplr(x_values_stimulus)], [mean_AC_stimulus_information - sem_AC_stimulus, ...
    fliplr(mean_AC_stimulus_information + sem_AC_stimulus)], 'r', 'FaceAlpha', 0.2, 'EdgeColor', 'none');
xline(7/samplingRate, '--', 'LineWidth', 1.5, 'Color', 'k');
hold off;
xlabel('Time (s)');
ylabel('Information (bits)');
title('Cumulative Stimulus Category Information');
legend([p3, p4], 'Location', 'Best');
ylim([0, 0.6]);

total_time_end = toc(total_time_start);
fprintf('Total Calculation Time =: %.4f Sekunden\n', total_time_end);

% Helper Function for mean and sem calculation 
function [mean_val, sem_val] = compute_stats(data)
    mean_val = mean(data, 1);
    std_val = std(data, 1);
    n = size(data, 1);
    sem_val = std_val / sqrt(n);
end

