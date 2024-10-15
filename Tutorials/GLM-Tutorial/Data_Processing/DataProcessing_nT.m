clc, clear, warning off,

%% Load data, split, normalize/z-Score 

load ('aligned_data_simple.mat');

PPC_data = struct();
AC_data = struct();
index_AC = 1;
index_PPC = 1;

for dataset = 1:length(aligned_data)
    name = aligned_data(dataset).dataset;
    if mod(dataset, 2) == 0
        
        neural_data_stimulus = aligned_data(dataset).r_aligned_stimulus;
        neural_data_choice = aligned_data(dataset).r_aligned_choice;
        [numTrials, numCells, ~] = size(neural_data_stimulus);

        % Normalize and Z-Score Data
        for cellIdx = 1:numCells
            for trialIdx = 1:numTrials
                maxVal_stimulus = max(neural_data_stimulus(trialIdx, cellIdx, :));
                maxVal_choice = max(neural_data_choice(trialIdx, cellIdx, :));

                if maxVal_stimulus ~= 0
                    neural_data_stimulus(trialIdx, cellIdx, :) = neural_data_stimulus(trialIdx, cellIdx, :) / maxVal_stimulus;
                end
                if maxVal_choice ~= 0
                    neural_data_choice(trialIdx, cellIdx, :) = neural_data_choice(trialIdx, cellIdx, :) / maxVal_choice;
                end
            end
        end
        neural_data_stimulus = zscore(neural_data_stimulus,[], 2);
        neural_data_choice = zscore(neural_data_choice,[], 2);
        PPC_data(index_PPC).data = name;
        PPC_data(index_PPC).neural_data_stimulus = neural_data_stimulus;
        PPC_data(index_PPC).neural_data_choice = neural_data_choice;
        PPC_data(index_PPC).x_pos_choice = aligned_data(dataset).x_pos_aligned_choice;
        PPC_data(index_PPC).y_pos_choice = aligned_data(dataset).y_pos_aligned_choice;
        PPC_data(index_PPC).x_pos_stimulus = aligned_data(dataset).x_pos_aligned_stimulus;
        PPC_data(index_PPC).y_pos_stimulus = aligned_data(dataset).y_pos_aligned_stimulus;
        PPC_data(index_PPC).choice = aligned_data(dataset).is_left_choice;
        PPC_data(index_PPC).sound_location = ceil(aligned_data(dataset).sound_location / 2);
        PPC_data(index_PPC).stimulus = aligned_data(dataset).is_left_stimulus;
        index_PPC = index_PPC+1;
    else
        neural_data_stimulus = aligned_data(dataset).r_aligned_stimulus;
        neural_data_choice = aligned_data(dataset).r_aligned_choice;
        [numTrials, numCells, ~] = size(neural_data_stimulus);

        % Normalize and Z-Score Data
        for cellIdx = 1:numCells
            for trialIdx = 1:numTrials
                maxVal_stimulus = max(neural_data_stimulus(trialIdx, cellIdx, :));
                maxVal_choice = max(neural_data_choice(trialIdx, cellIdx, :));

                if maxVal_stimulus ~= 0
                    neural_data_stimulus(trialIdx, cellIdx, :) = neural_data_stimulus(trialIdx, cellIdx, :) / maxVal_stimulus;
                end
                if maxVal_choice ~= 0
                    neural_data_choice(trialIdx, cellIdx, :) = neural_data_choice(trialIdx, cellIdx, :) / maxVal_choice;
                end
            end
        end
        neural_data_stimulus = zscore(neural_data_stimulus,[], 2);
        neural_data_choice = zscore(neural_data_choice,[], 2);
        AC_data(index_AC).data = name;
        AC_data(index_AC).neural_data_stimulus = neural_data_stimulus;
        AC_data(index_AC).neural_data_choice = neural_data_choice;
        AC_data(index_AC).x_pos_choice = aligned_data(dataset).x_pos_aligned_choice;
        AC_data(index_AC).y_pos_choice = aligned_data(dataset).y_pos_aligned_choice;
        AC_data(index_AC).x_pos_stimulus = aligned_data(dataset).x_pos_aligned_stimulus;
        AC_data(index_AC).y_pos_stimulus = aligned_data(dataset).y_pos_aligned_stimulus;
        AC_data(index_AC).choice = aligned_data(dataset).is_left_choice;
        AC_data(index_AC).sound_location = ceil(aligned_data(dataset).sound_location / 2);
        AC_data(index_AC).stimulus = aligned_data(dataset).is_left_stimulus;
        index_AC = index_AC+1;
    end
   
end

%% Split Training and Test Data AC_Data

k_split = 5;
trainIndices_AC = struct();
testIndices_AC = struct();

for data = 1:length(AC_data)
    name = AC_data(data).data;
    [numTrials, numCells, ~] = size(AC_data(data).neural_data_stimulus);
    Trials = 1:1:numTrials;
    cvPartition = cvpartition(Trials, 'KFold', k_split);

    for fold = 1:k_split
        field = ['Fold_', num2str(fold)];
        testIndices_AC(data).(field) =  find(cvPartition.test(fold));
    end
end


%% Split Training and Test Data PPC_Data

k_split = 5;
trainIndices_PPC = struct();
testIndices_PPC = struct();

for data = 1:length(PPC_data)
    name = PPC_data(data).data;
    [numTrials, numCells, ~] = size(PPC_data(data).neural_data_stimulus);
    Trials = 1:1:numTrials;
    cvPartition = cvpartition(Trials, 'KFold', k_split);

    for fold = 1:k_split
        field = ['Fold_', num2str(fold)];
        testIndices_PPC(data).(field) =  find(cvPartition.test(fold));
    end
end



for data_idx = 1:length(PPC_data)    
     PPC_data(data_idx).testIndices = testIndices_PPC(data_idx);
     AC_data(data_idx).testIndices = testIndices_AC(data_idx);
end


save('AC_data_nt.mat', 'AC_data');
save('PPC_data_nt.mat', 'PPC_data');
