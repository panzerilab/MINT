clc, clear,warning off,
rng(123);

% Load Data
load('AC_data.mat');
load('PPC_data.mat');

[~,nCc_AC,nTP_choice] = size(AC_data(6).neural_data_choice);
[~,nCs_AC,nTP_stimulus] = size(AC_data(6).neural_data_stimulus);
[~,nCs_PPC,~] = size(AC_data(7).neural_data_stimulus);

info_out_choice_AC = [];
info_out_choice_PPC = [];
info_out_stimulus_AC = [];

info_out_choice_AC_all = zeros(nCs_AC, nTP_choice);
info_out_choice_PPC_all = zeros(nCs_PPC,nTP_choice);
info_out_stimulus_AC_all = zeros(nCs_AC, nTP_stimulus);

opts.distribution = 'binomial';
opts.link = 'logit';
kFolds = 5;
MI_opts = struct('verbose', false, 'method', 'dr', 'bias', 'pt');

% -------------------------------------------------------------------------
% 1. Single Neuron Stimulus and Choice Information in AC
% -------------------------------------------------------------------------

disp('Calculating Single Neuron Stimulus and Choice Information in AC');

neural_data_choice = AC_data(6).neural_data_choice;
neural_data_stimulus = AC_data(6).neural_data_stimulus;
choice = AC_data(6).choice;
stimulus = AC_data(6).stimulus;
[numTrials, numCells, numTimePoints_choice] = size(neural_data_choice);
[~,~,numTimePoints_stimulus] = size(neural_data_stimulus);
predicted_choice_AC = zeros(numTimePoints_choice, numTrials);
predicted_stimulus_AC = zeros(numTimePoints_stimulus, numTrials);

for cellIdx = 1:numCells
    disp(['Neuron: ', num2str(cellIdx)]);
    neuron_idx = sprintf('Neuron_%d', cellIdx);

    for tPc = 1:numTimePoints_choice
        for i = 1:length(fields(AC_data(6).trainIndices))
            fold_idx = ['Fold_', num2str(i)];
            test_idxs = AC_data(6).testIndices.(fold_idx);
            neuralDataAtTime = neural_data_choice(:, cellIdx, tPc);
            [tmpPredLabels] = fitglm_pipeline(neuralDataAtTime, choice, test_idxs, opts);
            predicted_choice_AC(tPc, test_idxs)= tmpPredLabels;
        end
    end
    predicted_choice_AC = predicted_choice_AC';

    for tPs = 1:numTimePoints_stimulus
        for i = 1:length(fields(AC_data(6).trainIndices))
            fold_idx = ['Fold_', num2str(i)];
            test_idxs = AC_data(6).testIndices.(fold_idx);
            neuralDataAtTime = neural_data_stimulus(:, cellIdx, tPs);
            [tmpPredLabels] = fitglm_pipeline(neuralDataAtTime, stimulus, test_idxs, opts);
            predicted_stimulus_AC(tPs, test_idxs)= tmpPredLabels;
        end
    end
    predicted_stimulus_AC = predicted_stimulus_AC';


    for tPc = 1:numTimePoints_choice
        predicted_choice_tmp = double(mean(predicted_choice_AC(:, 1:tPc), 2) >= 0.5);
        info_out = information(predicted_choice_tmp', choice', MI_opts, {'I', 'Ish'});
        info_out_choice_AC(tPc) =  info_out{1, 1};
    end
    info_out_choice_AC_all(cellIdx,:) = info_out_choice_AC;

    for tPs = 1:numTimePoints_stimulus
        predicted_stimulus_tmp = double(mean(predicted_stimulus_AC(:, 1:tPs), 2) >= 0.5);
        info_out = information(predicted_stimulus_tmp', stimulus', MI_opts, {'I', 'Ish'});
        info_out_stimulus_AC(tPs) =  info_out{1, 1};
    end
    info_out_stimulus_AC_all(cellIdx,:) = info_out_stimulus_AC;
end


% mean_AC_choice_information = mean(info_out_choice_AC_all, 1);
% std_AC_choice_information = std(info_out_choice_AC_all,1);
% mean_AC_stimulus_information = mean(info_out_stimulus_AC_all, 1);
% std_AC_stimulus_information = std(info_out_stimulus_AC_all,1);

% -------------------------------------------------------------------------
% 2. Single Neuron Choice Information in PPC
% -------------------------------------------------------------------------

disp('Calculating Single Neuron Choice Information in PPC');
neural_data_choice = PPC_data(7).neural_data_choice;
choice = PPC_data(7).choice;
[numTrials, numCells, numTimePoints_choice] = size(neural_data_choice);
predicted_choice_AC = zeros(numTimePoints_choice, numTrials);

for cellIdx = 1:numCells
    disp(['Neuron: ', num2str(cellIdx)]);
    neuron_idx = sprintf('Neuron_%d', cellIdx);

    for tPc = 1:numTimePoints_choice
        for i = 1:length(fields(PPC_data(7).trainIndices))
            fold_idx = ['Fold_', num2str(i)];
            test_idxs = PPC_data(7).testIndices.(fold_idx);
            neuralDataAtTime = neural_data_choice(:, cellIdx, tPc);
            [tmpPredLabels] = fitglm_pipeline(neuralDataAtTime, choice, test_idxs, opts);
            predicted_choice_PPC(tPc, test_idxs)= tmpPredLabels;
        end
    end
    predicted_choice_PPC = predicted_choice_PPC';

    for tPc = 1:numTimePoints_choice
        predicted_choice_tmp = double(mean(predicted_choice_PPC(:, 1:tPc), 2) >= 0.5);
        info_out = information(predicted_choice_tmp', choice', MI_opts, {'I', 'Ish'});
        info_out_choice_PPC(tPc) =  info_out{1, 1};
    end
    info_out_choice_PPC_all(cellIdx,:) = info_out_choice_PPC;
end

% mean_PPC_choice_information = mean(info_out_choice_PPC_all, 1);
% std_PPC_choice_information = std(info_out_choice_PPC_all,1);

save('A1_B77_0912_choice_singleNeuronInformation.mat', 'info_out_choice_AC_all');
save('A1_B77_0912_stimulus_singleNeuronInformation.mat', 'info_out_stimulus_AC_all');
save('PPC_B67_0629_choice_singleNeuronInformation.mat', 'info_out_choice_PPC_all');

