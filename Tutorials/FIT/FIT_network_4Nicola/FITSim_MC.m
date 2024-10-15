clear all;
rng(0);

tic;
simReps = 10; % number of repetitions of the simulation

nPopulations = 4;
nSubpopulations = 2;
nTrials = 200;
nTimepoints = 30;
delay = 5;
stim_window = [3 12];
stim_window_2 = [3+delay+2 12+delay+2];
nShuff = 20; % number of shufflings
null_samples = 200; % number of null hypothesis samples obtained by bootstrapping shufflings across repetitions

c = 1;

X = zeros(nPopulations, nSubpopulations, nTimepoints, nTrials);
infoS1 = zeros(nPopulations,nTimepoints, simReps);
infoS2 = zeros(nPopulations,nTimepoints, simReps);

FIT1 = zeros(nPopulations,nPopulations, simReps);
FIT2 = zeros(nPopulations,nPopulations, simReps);

delay_shifts = [0:1:nTimepoints-1];
FIT1_shift = nan(nTimepoints, length(delay_shifts), simReps);
FIT2_shift = nan(nTimepoints, length(delay_shifts), simReps);

FIT1_sh.conditioned = zeros(nPopulations,nPopulations, simReps, nShuff);
FIT2_sh.conditioned = zeros(nPopulations,nPopulations, simReps, nShuff);
TE_sh.conditioned = zeros(nPopulations,nPopulations, simReps, nShuff);
FIT1_sh.simple = zeros(nPopulations,nPopulations, simReps, nShuff);
FIT2_sh.simple = zeros(nPopulations,nPopulations, simReps, nShuff);
TE_sh.simple = zeros(nPopulations,nPopulations, simReps, nShuff);

TE = zeros(nPopulations,nPopulations, simReps);


opts = [];
opts.verbose = false;
opts.method = "dr";
opts.bias = 'qe';
opts.xtrp = 10;
opts.btsp = 0;
opts.n_binsS = 2; % Number of stimulus values
opts.n_binsX = 3;
opts.n_binsY = 3;
opts.bin_methodS = 'none';
opts.bin_methodX = 'eqpop';
opts.taux = -delay;
opts.tauy = -delay;
opts.btsp_variables = {};
opts.parallel = 0;
opts.singleTimepoint = true;
for repIdx = 1:simReps
    repIdx
    X = zeros(nPopulations, nSubpopulations, nTimepoints, nTrials);
    S1 = randi(opts.n_binsS,1,nTrials);
    S2 = randi(opts.n_binsS,1,nTrials);

    epsX = 0.5;   

    for t=1:nTimepoints

        X(:,:,t,:) = epsX*randn(nPopulations, nSubpopulations,1, nTrials);
      
        if (t > stim_window(1))  && (t <  stim_window(2))
            X(1,1,t,:) = squeeze(X(1,1,t,:)) +epsX*randn+2*(S1-1.5)';
            X(2,2,t,:) = squeeze(X(2,2,t,:)) +epsX*randn+2*(S2-1.5)';
        end

        if (t > stim_window_2(1))  && (t <  stim_window_2(2))
            X(4,1,t,:) = squeeze(X(4,1,t,:)) +epsX*randn+2*(S1-1.5)';
        end
 
        if t > delay 
            X(1,2,t,:) = X(1,2,t,:) + c*X(2,2,t-delay,:);
            X(2,1,t,:) = X(2,1,t,:) + c*X(1,1,t-delay,:);
            X(3,1,t,:) = X(3,1,t,:) + c*X(1,1,t-delay,:);
            X(4,2,t,:) = X(4,2,t,:) + 3*c*X(3,2,t-delay,:); % noise transfer
        end
        
        opts.bin_methodY = 'none';

        for neuron = 1:nPopulations
            infoS1(neuron, t, repIdx) = cell2mat(information(squeeze(sum(X(neuron,:,t,:),2))', S1 ,opts,{'I'}));
            infoS2(neuron, t, repIdx) = cell2mat(information(squeeze(sum(X(neuron,:,t,:),2))', S2 ,opts,{'I'}));
        end

    end
    X = squeeze(sum(X,2));
    opts.bin_methodY = 'eqpop';
    % for d = 1:length(delay_shifts)
    %     for t=1:nTimepoints
    %         delay_shift = delay_shifts(d);
    %         if t-delay_shift > 0
    %             X_pres = squeeze(X(:,t,:));
    %             X_past = squeeze(X(:,t-delay_shift,:));
    %             FIT1_shift(t, d, repIdx) = FIT(S1,X_past(1,:),X_past(2,:),X_pres(2,:),opts);
    %             FIT2_shift(t, d, repIdx) = FIT(S2,X_past(2,:),X_past(1,:),X_pres(1,:),opts);
    %         end
    %     end
    % end 

    t = 7+delay;

    X_pres = squeeze(X(:,t,:)); 
    X_past = squeeze(X(:,t-delay,:));

    for neuron1 = 1:nPopulations
       for neuron2 = 1:nPopulations
           if neuron1 ~= neuron2
               FIT1(neuron1, neuron2, repIdx) = FIT(S1,X_past(neuron1,:),X_past(neuron2,:),X_pres(neuron2,:),opts);            
               FIT2(neuron1, neuron2, repIdx) = FIT(S2,X_past(neuron1,:),X_past(neuron2,:),X_pres(neuron2,:),opts);                         
               TE(neuron1, neuron2, repIdx) = cell2mat(transferentropy2(X(neuron1,1:t,:), X(neuron2,1:t,:),opts, {'TE'}));          
           end 
       end 
    end 
 
    % Compute null hypothesis for FIT and TE
    for shIdx = 1:nShuff
        shIdx
        shuffledS1 = S1(randperm(numel(S1)));
        shuffledS2 = S2(randperm(numel(S2)));
        S1_S2 = map_Nd_array_to_1d([S1;S2]);

        for neuron1 = 1:nPopulations

            % Permute data to have trials in the first dimension (as
            % required by the shuffle function)
            sender = X(neuron1,:,:);
            sender_perm = permute(sender,[3 1 2]);

            % constrained shuffling of fixed values of S1 and S2
            sender_s1_s2_shuff_perm = shuffle(S1_S2',sender_perm,1,[1 0 0 ]);
  %          sender_s1_s2_shuff_2 = shuffle(S1_S2',sender,1,[0 0 1 ]); TBD
            sender_s1_s2_shuff = permute(sender_s1_s2_shuff_perm,[2 3 1]);

            % Unconstrained shuffling of the sender
            sender_all_shuff_perm = shuffle(S1',sender_perm,0,[1 0 0 ]);
            sender_all_shuff = permute(sender_all_shuff_perm,[2 3 1]);

            % Shuffle the sender at fixed values of S1
            sender_s1_shuff_perm = shuffle(S1',sender_perm,1,[1 0 0 ]);
            sender_s1_shuff = permute(sender_s1_shuff_perm,[2 3 1]);

            % Shuffle the sender at fixed values of S2
            sender_s2_shuff_perm = shuffle(S2',sender_perm,1,[1 0 0 ]);
            sender_s2_shuff = permute(sender_s2_shuff_perm,[2 3 1]);

            X_s1_shuff = X; 
            X_s2_shuff = X; 
            X_all_shuff = X;
            X_s1_s2_shuff = X; 
            X_s1_shuff(neuron1,:,:) = sender_s1_shuff;
            X_s2_shuff(neuron1,:,:) = sender_s2_shuff;
            X_s1_s2_shuff(neuron1,:,:) = sender_s1_s2_shuff;
            X_all_shuff(neuron1,:,:) = sender_all_shuff;

            % compute past and present of neural signals
            X_all_shuff_pres = squeeze(X_all_shuff(:,t,:)); 
            X_all_shuff_past = squeeze(X_all_shuff(:,t-delay,:)); 
            X_s1_shuff_pres = squeeze(X_s1_shuff(:,t,:)); 
            X_s1_shuff_past = squeeze(X_s1_shuff(:,t-delay,:)); 
            X_s2_shuff_pres = squeeze(X_s2_shuff(:,t,:)); 
            X_s2_shuff_past = squeeze(X_s2_shuff(:,t-delay,:)); 

            for neuron2 = 1:nPopulations
                if neuron1 ~= neuron2
                    % Conditioned shuffling (fixed stimulus)FIT(S1,X_past(neuron1,:),X_past(neuron2,:),X_pres(neuron2,:),opts);   
                    FIT1_sh.conditioned(neuron1, neuron2, repIdx, shIdx) = FIT(S1,X_s1_shuff_past(neuron1,:),X_s1_shuff_past(neuron2,:),X_s1_shuff_pres(neuron2,:),opts);            
                    FIT2_sh.conditioned(neuron1, neuron2, repIdx, shIdx) = FIT(S2,X_s2_shuff_past(neuron1,:),X_s2_shuff_past(neuron2,:),X_s2_shuff_pres(neuron2,:),opts);
                    TE_sh.conditioned(neuron1, neuron2, repIdx, shIdx) = cell2mat(transferentropy2(X_s1_s2_shuff(neuron1,1:t,:),X_s1_s2_shuff(neuron2,1:t,:),opts, {'TE'}));
                    
                    % Simple S shuffling
                    FIT1_sh.simple(neuron1, neuron2, repIdx, shIdx) = FIT(shuffledS1,X_past(neuron1,:),X_past(neuron2,:),X_pres(neuron2,:),opts);            
                    FIT2_sh.simple(neuron1, neuron2, repIdx, shIdx) = FIT(shuffledS2,X_past(neuron1,:),X_past(neuron2,:),X_pres(neuron2,:),opts);
                    TE_sh.simple(neuron1, neuron2, repIdx, shIdx) = cell2mat(transferentropy2(X_all_shuff(neuron1,1:t,:), X_all_shuff(neuron2,1:t,:),opts, {'TE'}));    % change the shuffling                        
                end 
            end 
        end
    end 
end

infoS1_t = zeros(4,simReps);
infoS2_t = zeros(4,simReps);

for neuron = 1:nPopulations
    if neuron == 1 || neuron == 4
        t = 5;
    elseif neuron == 2|| neuron == 3 || neuron == 5
        t = 12;
    end 
    infoS1_t(neuron,:) = squeeze(infoS1(neuron, t, :));
end 

for neuron = 1:nPopulations
    if neuron == 2 
        t = 5;
    else 
        t = 12;
    end 
    infoS2_t(neuron,:) = squeeze(infoS2(neuron, t, :));
end 

FIT1_mean = mean(FIT1, 3);
FIT2_mean = mean(FIT2, 3);
TE_mean = mean(TE, 3);
%%

% Bootstrap shuffled values across repetitions of the simulations to boost
% the number of null hypothesis samples
FIT1_sh_perm.simple = permute(FIT1_sh.simple,[3 1 2 4]);
FIT1_sh_mean.simple = btstrp_shuff(FIT1_sh_perm.simple,null_samples);
FIT2_sh_perm.simple = permute(FIT2_sh.simple,[3 1 2 4]);
FIT2_sh_mean.simple = btstrp_shuff(FIT2_sh_perm.simple,null_samples);

FIT1_sh_perm.conditioned = permute(FIT1_sh.conditioned,[3 1 2 4]);
FIT1_sh_mean.conditioned = btstrp_shuff(FIT1_sh_perm.conditioned,null_samples);
FIT2_sh_perm.conditioned = permute(FIT2_sh.conditioned,[3 1 2 4]);
FIT2_sh_mean.conditioned = btstrp_shuff(FIT2_sh_perm.conditioned,null_samples);

% Take the element-wise maximum across the two null hypotheses
FIT1_sh_mean.max = max(FIT1_sh_mean.simple,FIT1_sh_mean.conditioned);
FIT2_sh_mean.max = max(FIT2_sh_mean.simple,FIT2_sh_mean.conditioned);


TE_sh_perm.simple = permute(TE_sh.simple,[3 1 2 4]);
TE_sh_mean.simple = btstrp_shuff(TE_sh_perm.simple,null_samples);

TE_sh_perm.conditioned = permute(TE_sh.conditioned,[3 1 2 4]);
TE_sh_mean.conditioned = btstrp_shuff(TE_sh_perm.conditioned,null_samples);

TE_sh_mean.max = max(TE_sh_mean.simple,TE_sh_mean.conditioned);

% TE_sh_perm = permute(TE_sh,[3 1 2 4]);
% TE_sh_mean = btstrp_shuff(TE_sh_perm,100);

% Determine significance thresholds for FIT 1 and FIT 2
sig_thresh_FIT1 = prctile(FIT1_sh_mean.max,99,3);
sig_thresh_FIT2 = prctile(FIT2_sh_mean.max,99,3);
sig_thresh_TE = prctile(TE_sh_mean.max,99,3);

toc

%save('FIT_network_results_2006')
%%
for neuron = 1:nPopulations

    figure;
    meanS1 = mean(infoS1_t(neuron, :));
    meanS2 = mean(infoS2_t(neuron, :));
    % semS1_curr = semS1(neuron);
    % semS2_curr = semS2(neuron);

    colorS1 = [123 148 146]/255;
    colorS2 = [190 161 127]/255;

    bar(1, meanS1, 'FaceColor', colorS1); 
    hold on;
    bar(2, meanS2, 'FaceColor', colorS2);
    % errorbar(1, meanS1, semS1_curr, 'k', 'linestyle', 'none');
    % errorbar(2, meanS2, semS2_curr, 'k', 'linestyle', 'none');

    set(gca, 'XTick', [1 2], 'XTickLabel', {'', ''});
    yLimits = [0, 0.65]; 
    ylim(yLimits);

    set(gca, 'YColor', 'none');
    set(gca, 'YTickLabel', []);
    set(gca, 'Box', 'off');

    %title(['Neuron ' num2str(neuron)]);
    saveas(gcf, ['Information_FITsim_Neuron' num2str(neuron) '.svg']);
    close(gcf);
end
%%

infoS1_mean_1 = squeeze(mean(infoS1(1,:,:),3));
% figure('Position',[270,197,798,277])
% subplot(1,2,1)
plot(infoS1_mean_1,'LineWidth',2)
xlabel('time')
ylabel('[bits]')
title('info S1')

infoS1_mean_2 = squeeze(mean(infoS1(2,:,:),3));
% subplot(1,2,2)
plot(infoS1_mean_1,'LineWidth',2)
xlabel('time')
ylabel('[bits]')
title('info S1')

% Plot TE and FIT
figure('Position',[25,230,1179,296])
subplot(1,3,1)
imagesc(TE_mean)
for i = 1:nPopulations
    for j = 1:nPopulations
        if(TE_mean(i,j)>sig_thresh_TE(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
yticks(1:nPopulations)
yticklabels_array = cellstr(num2str((1:nPopulations)', '%d'));
title('TE')
colorbar()

subplot(1,3,2)
imagesc(FIT1_mean)
for i = 1:nPopulations
    for j = 1:nPopulations
        if(FIT1_mean(i,j)>sig_thresh_FIT1(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
title('FIT_{S1}')
yticks(1:nPopulations)
yticklabels_array = cellstr(num2str((1:nPopulations)', '%d'));
colorbar()

subplot(1,3,3)
imagesc(FIT2_mean)
for i = 1:nPopulations
    for j = 1:nPopulations
        if(FIT2_mean(i,j)>sig_thresh_FIT2(i,j))
            text(j,i,'*')
        end
    end
end
xlabel('receiver')
ylabel('sender')
yticks(1:nPopulations)
yticklabels_array = cellstr(num2str((1:nPopulations)', '%d'));
title('FIT_{S2}')
colorbar()

saveas(gcf, 'TE_and_FIT_FITsim.png')

%%
% figure;
% 
% % Datenvorbereitung für FIT1
% infoS1_mean_1 = squeeze(mean(infoS1(1,:,:), 3));
% infoS1_mean_2 = squeeze(mean(infoS1(2,:,:), 3));
% SEM_1 = squeeze(std(infoS1(1,:,:), 0, 3) / sqrt(size(infoS1, 3)));
% SEM_2 = squeeze(std(infoS1(2,:,:), 0, 3) / sqrt(size(infoS1, 3)));
% FIT1_Shift_mean = mean(FIT1_shift,3);
% y_limits_1 = [min(min(infoS1_mean_1 - SEM_1), min(infoS1_mean_2 - SEM_2)) max(max(infoS1_mean_1 + SEM_1), max(infoS1_mean_2 + SEM_2))+0.1];
% 
% % Datenvorbereitung für FIT2
% infoS2_mean_1 = squeeze(mean(infoS2(1,:,:), 3));
% infoS2_mean_2 = squeeze(mean(infoS2(2,:,:), 3));
% SEM_1_2 = squeeze(std(infoS2(1,:,:), 0, 3) / sqrt(size(infoS2, 3)));
% SEM_2_2 = squeeze(std(infoS2(2,:,:), 0, 3) / sqrt(size(infoS2, 3)));
% FIT2_Shift_mean = mean(FIT2_shift,3);
% y_limits_2 = [min(min(infoS2_mean_1 - SEM_1_2), min(infoS2_mean_2 - SEM_2_2)) max(max(infoS2_mean_1 + SEM_1_2), max(infoS2_mean_2 + SEM_2_2))+0.1];
% 
% % Gemeinsame y-Limits festlegen
% y_limits = [min(y_limits_1(1), y_limits_2(1)), max(y_limits_1(2), y_limits_2(2))];
% 
% % Erstellen des tiledLayouts
% t = tiledlayout(6, 2, 'TileSpacing', 'Compact', 'Padding', 'Compact');
% 
% % Plot 1 - FIT1, MI(X1;S1)
% nexttile
% colorS1 = [123 148 146]/255;
% plot(infoS1_mean_1, 'LineWidth', 1.5, 'Color', colorS1)
% hold on
% fill([1:length(infoS1_mean_1), length(infoS1_mean_1):-1:1], ...
%     [infoS1_mean_1 - SEM_1, fliplr(infoS1_mean_1 + SEM_1)], colorS1, ...
%     'FaceAlpha', 0.3, 'EdgeColor', 'none')
% hold off
% ylabel('[bits]', 'FontSize', 12, 'FontName', 'Arial')
% set(gca, 'YLim', y_limits, 'FontSize', 10, 'FontName', 'Arial')
% title('MI(X1;S1)')
% xticks([1, 5:5:25, 30])
% set(gca, 'XTickLabel', [], 'TickDir', 'out', 'Box', 'off')
% 
% % Plot 2 - FIT2, MI(X2;S2)
% nexttile
% colorS2 = [190 161 127]/255;
% plot(infoS2_mean_2, 'LineWidth', 1.5, 'Color', colorS2)
% hold on
% fill([1:length(infoS2_mean_2), length(infoS2_mean_2):-1:1], ...
%     [infoS2_mean_2 - SEM_2_2, fliplr(infoS2_mean_2 + SEM_2_2)], colorS2, ...
%     'FaceAlpha', 0.3, 'EdgeColor', 'none')
% hold off
% set(gca, 'YLim', y_limits,'YTick', [], 'TickDir', 'out', 'Box', 'off')
% title('MI(X2;S2)')
% xticks([1, 5:5:25, 30])
% set(gca, 'XTickLabel', [], 'TickDir', 'out', 'Box', 'off')
% 
% % Heatmap - FIT1
% nexttile([4 1])
% h = imagesc(1:nTimepoints, delay_shifts, FIT1_Shift_mean');
% alpha_data = ~isnan(FIT1_Shift_mean');
% set(h, 'AlphaData', alpha_data);
% ylabel(' Δt', 'FontSize', 12, 'FontName', 'Arial')
% colormap(jet);
% yticks([0, 5:5:25, 29])
% xticks([1, 5:5:25, 30])
% set(gca, 'XTickLabel', [], 'FontSize', 10, 'FontName', 'Arial')
% 
% % Heatmap - FIT2
% nexttile([4 1])
% h = imagesc(1:nTimepoints, delay_shifts, FIT2_Shift_mean');
% alpha_data = ~isnan(FIT2_Shift_mean');
% set(h, 'AlphaData', alpha_data);
% colorbar
% set(get(colorbar, 'Title'), 'String', '[bits]', 'FontSize', 12, 'FontName', 'Arial');
% colormap(jet);
% yticks([0, 5:5:25, 29])
% xticks([1, 5:5:25, 30])
% set(gca, 'XTickLabel', [], 'YTickLabel', [], 'FontSize', 10, 'FontName', 'Arial')
% 
% 
% % Plot 3 - FIT1, MI(X2;S1)
% nexttile
% plot(infoS1_mean_2, 'LineWidth', 1.5, 'Color', colorS1)
% hold on
% fill([1:length(infoS1_mean_2), length(infoS1_mean_2):-1:1], ...
%     [infoS1_mean_2 - SEM_2, fliplr(infoS1_mean_2 + SEM_2)], colorS1, ...
%     'FaceAlpha', 0.3, 'EdgeColor', 'none')
% hold off
% xlabel('t', 'FontSize', 12, 'FontName', 'Arial')
% ylabel('[bits]', 'FontSize', 12, 'FontName', 'Arial')
% set(gca, 'YLim', y_limits, 'FontSize', 10, 'FontName', 'Arial')
% title('MI(X2;S1)')
% xticks([1, 5:5:25, 30])
% set(gca, 'TickDir', 'out', 'Box', 'off')
% 
% % Plot 4 - FIT2, MI(X1;S2)
% nexttile
% plot(infoS2_mean_1, 'LineWidth', 1.5,'Color', colorS2)
% hold on
% fill([1:length(infoS2_mean_1), length(infoS2_mean_1):-1:1], ...
%     [infoS2_mean_1 - SEM_1_2, fliplr(infoS2_mean_1 + SEM_1_2)], colorS2, ...
%     'FaceAlpha', 0.3, 'EdgeColor', 'none')
% hold off
% xlabel('t', 'FontSize', 12, 'FontName', 'Arial')
% set(gca, 'YLim', y_limits,'YTick', [], 'TickDir', 'out', 'Box', 'off')
% title('MI(X1;S2)')
% xticks([1, 5:5:25, 30])
% set(gca, 'TickDir', 'out', 'Box', 'off')
% linkaxes([t.Children], 'x');
% 
% saveas(gcf, 'Combined_DelaySweep_FIT1_FIT2.svg');
% 
