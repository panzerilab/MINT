clc, clear,warning off,
rng(123);

MI_opts.verbose = false;
MI_opts.method = "dr";
MI_opts.bias = 'naive';
MI_opts.n_binsX = 3;
MI_opts.btsp = 100;
MI_opts.bin_methodX = 'eqpop';
MI_opts.bin_methodY = 'none';
% 
II_opts.bias = 'naive';
II_opts.bin_methodR = 'eqpop';
II_opts.n_binsR = 3;
II_opts.btsp = 100;


 %
% AC_data_training = load('AC_data_uncut.mat');
% PPC_data_training = load('PPC_data_uncut.mat');
AC_data = load('AC_data_unbinned_cut.mat');
PPC_data = load('PPC_data_unbinned_cut.mat');
AC_data = AC_data.AC_data;
PPC_data = PPC_data.PPC_data;

%% AC Information

datasets_AC = {'A1_B67_0630', 'A1_B67_0703', 'A1_B68_0702', 'A1_B75_0826', 'A1_B75_0908', 'A1_B77_0905', 'A1_B77_0912'};
information_stim_AC = [];
information_choice_AC = [];
information_II_AC = [];
MI_breakdown_AC_stim = struct();
MI_breakdown_AC_choice = struct();
MI_breakdown_AC_II = struct();


for dataset = 1:length(datasets_AC)
    datasetName = datasets_AC{dataset};
    disp([num2str(dataset),'.Dataset: ', datasetName]);
    data = AC_data.(datasetName);
    numCells = size(data(1).response,1);
    numTrials = length(data);
    R_tot = [];
    S = [];
    C = [];

    for i = 1:numTrials
        data(i).stimulus = data(i).sound_loc;
        % Calculate mean response each trial 
        for cell = 1:numCells
            
            data(i).mean_response = mean(data(i).response,2);
        end
    end

    for i = 1:numTrials
        R_tot = [R_tot, data(i).mean_response];
        S = [S, data(i).stimulus];
        C = [C, data(i).left_turn];
    end

    pairIdx = 0;
    for cell = 1:numCells
        infoBdw_stim = information(R_tot(cell,:),S,MI_opts,{'I'});
        information_out_stim = infoBdw_stim{1, 1}(1)-(mean(infoBdw_stim{1, 1}(2:end)));
        MI_breakdown_AC_stim.(datasetName)(cell) = information_out_stim;

        infoBdw_choice = information(R_tot(cell,:),C,MI_opts,{'I'});
        information_out_choice = infoBdw_choice{1, 1}(1)-(mean(infoBdw_choice{1, 1}(2:end)));
        MI_breakdown_AC_choice.(datasetName)(cell) = information_out_choice; 
        
        information_out_II = Intersection_Information(S, R_tot(cell,:), C, II_opts);
        information_out_II_value = information_out_II.biased.value;
        if information_out_II_value < 0 
            information_out_II_value = 0;
        end 
        MI_breakdown_AC_II.(datasetName)(cell) = information_out_II_value;
    end

    information_stim_AC = [information_stim_AC, mean(MI_breakdown_AC_stim.(datasetName))];
    information_choice_AC = [information_choice_AC, mean(MI_breakdown_AC_choice.(datasetName))];
    information_II_AC = [information_II_AC, mean(MI_breakdown_AC_II.(datasetName)')];
end

%% PPC Stimulus Information 

datasets_PPC = {'PPC_B67_0629', 'PPC_B67_0710', 'PPC_B67_0716', 'PPC_B68_0628', 'PPC_B68_0630', 'PPC_B68_0703', 'PPC_B77_0828'};
information_stim_PPC = [];
information_choice_PPC = [];
information_II_PPC = [];
MI_breakdown_PPC_stim = struct();
MI_breakdown_PPC_choice = struct();
MI_breakdown_PPC_II = struct();

for dataset = 1:length(datasets_PPC)
    datasetName = datasets_PPC{dataset};
    disp([num2str(dataset),'.Dataset: ', datasetName]);
    data = PPC_data.(datasetName);
    numCells = size(data(1).response,1);
    numTrials = length(data);
    R_tot = [];
    S = [];
    C = [];

    for i = 1:numTrials
        data(i).stimulus = data(i).sound_loc;    
        % Calculate mean response each trial 
        for cell = 1:numCells
            data(i).mean_response = mean(data(i).response,2);
        end
    end

    for i = 1:numTrials
        R_tot = [R_tot, data(i).mean_response];
        S = [S, data(i).stimulus];
        C = [C, data(i).left_turn];
    end


    pairIdx = 0;
    for cell = 1:numCells
        infoBdw_stim = information(R_tot(cell,:),S,MI_opts,{'I'});
        information_out_stim = infoBdw_stim{1, 1}(1)-(mean(infoBdw_stim{1, 1}(2:end)));
        MI_breakdown_PPC_stim.(datasetName)(cell) = information_out_stim;

        infoBdw_choice = information(R_tot(cell,:),C,MI_opts,{'I'});
        information_out_choice = infoBdw_choice{1, 1}(1)-(mean(infoBdw_choice{1, 1}(2:end)));
        MI_breakdown_PPC_choice.(datasetName)(cell) = information_out_choice;
       
       information_out_II = Intersection_Information(S, R_tot(cell,:), C, II_opts);
        information_out_II_value = information_out_II.biased.value;
        if information_out_II_value < 0 
            information_out_II_value = 0;
        end 
        MI_breakdown_PPC_II.(datasetName)(cell) =  information_out_II_value;
    end

    information_stim_PPC = [information_stim_PPC, mean(MI_breakdown_PPC_stim.(datasetName))];
    information_choice_PPC = [information_choice_PPC, mean(MI_breakdown_PPC_choice.(datasetName))];
    information_II_PPC = [information_II_PPC, mean(MI_breakdown_PPC_II.(datasetName))];
end



%% Plot
mean_choice_AC = mean(information_choice_AC);
sem_choice_AC = std(information_choice_AC)/sqrt(length(information_choice_AC));

mean_choice_PPC = mean(information_choice_PPC);
sem_choice_PPC = std(information_choice_PPC)/sqrt(length(information_choice_PPC));

mean_stim_AC = mean(information_stim_AC);
sem_stim_AC = std(information_stim_AC)/sqrt(length(information_stim_AC));

mean_stim_PPC = mean(information_stim_PPC);
sem_stim_PPC = std(information_stim_PPC)/sqrt(length(information_stim_PPC));

mean_II_AC = mean(information_II_AC);
sem_II_AC = std(information_II_AC)/sqrt(length(information_II_AC));

mean_II_PPC = mean(information_II_PPC);
sem_II_PPC = std(information_II_PPC)/sqrt(length(information_II_PPC));

data = [mean_stim_AC mean_stim_PPC; mean_choice_AC mean_choice_PPC; mean_II_AC mean_II_PPC];
errors = [sem_stim_AC sem_stim_PPC; sem_choice_AC sem_choice_PPC; sem_II_AC sem_II_PPC];

maxStim = max(mean_stim_AC + sem_stim_AC, mean_stim_PPC + sem_stim_PPC);
maxChoice = max(mean_choice_AC + sem_choice_AC, mean_choice_PPC + sem_choice_PPC);
maxII = max(mean_II_AC + sem_II_AC, mean_II_PPC + sem_II_PPC);

maxY = maxStim+0.005;

[~, p_Stim] = ttest(information_stim_AC,information_stim_PPC);
[~, p_Choice] = ttest(information_choice_AC, information_choice_PPC);
[~, p_II] = ttest(information_II_AC, information_II_PPC);

figure;
set(gcf, 'Position', [60, 60, 480, 240]); 

b = bar(data, 'grouped');
hold on;


b(1).FaceColor = '#cb181d'; 
b(2).FaceColor = '#2270b5'; 

b(1).EdgeColor = 'black';
b(2).EdgeColor = 'black';

b(1).LineWidth = 1.5; 
b(2).LineWidth = 1.5;

[numGroups, numBars] = size(data);
groupWidth = min(0.8, numBars/(numBars + 1.5));
for i = 1:numBars
    x = (1:numGroups) - groupWidth/2 + (2*i-1) * groupWidth / (2*numBars);
    e = errorbar(x, data(:,i), errors(:,i), 'k', 'linestyle', 'none');
    e.LineWidth = 1.5; 
end

ax = gca;
ax.FontSize = 14; 
ylim([0, maxY+0.01])
set(gca, 'XTick', 1:numGroups, 'XTickLabel', {'Stimulus', 'Choice', 'Intersection'});
ylabel('Information (bits)');
add_pvalue([p_Stim, p_Choice, p_II], maxY, [0.75,1.75,2.75], [1.25,2.25,3.25], 0.0025, [maxStim-maxY + 0.001, maxY-maxChoice - 0.01, maxY-maxII - 0.01])
hold off;

saveas(gcf, 'Information_plot_2_uncut.png');


%% Helper Function

function [mean_val, sem_val] = compute_stats(data)
    mean_val = mean(data);
    std_val = std(data);
    n = length(data);
    sem_val = std_val / sqrt(n);
end

function add_pvalue(p_values, y, x1, x2, text_heights, line_heights)
    for i = 1:numel(p_values)
        stars = get_stars(p_values(i));
        line([x1(i), x2(i)], [y - line_heights(i), y - line_heights(i)], 'Color', 'black', 'LineWidth', 1);
        text((x1(i) + x2(i)) / 2, y - line_heights(i) + text_heights, stars, 'HorizontalAlignment', 'center', 'Color', 'black');
    end
end

function stars = get_stars(p_value)
    if p_value < 0.001
        stars = '***';
    elseif p_value < 0.01
        stars = '**';
    elseif p_value < 0.05
        stars = '*';
    else
        stars = 'n.s.';
    end
end