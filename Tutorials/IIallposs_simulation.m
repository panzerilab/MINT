% Simulation to reproduce Fig. 3A results
% A population of 2 neurons encodes a bianry stimulus. We compute II and plot 
% stimulus and choice decision baoundaries for two fixed encoding/decoding angles.
% decoding an
clear all;

rng(0)

N_cells = 2;
N_trials_list = [10, 50, 100, 500, 1000]; %[10:10:100];
nfolds = 10;
repetitions = 5;
angle_range = [20,70]; %(0:2.5:90);
angle_range = angle_range*(pi/180);



info_amount = {'high_info'}; % info scenario name, not really needed here but kept for consistency with population simulation
mean_FR.high_info = 4;
delta_FR.high_info = 1;

% Info options
opts.bin_methodS = 'none';
opts.bin_methodR = 'none';
opts.bin_methodC = 'none';
opts.bin_methodX = 'none';
opts.bin_methodY = 'none';
opts.method = 'dr';
opts.bias = 'qe';
opts.xtrp = 10;
opts.verbose = 0;
opts_singlecell = opts;
opts_singlecell.bin_methodR = 'eqpop';
opts_singlecell.n_binsR = 2;

% SVM options
SVM_opts.optimize_params = false; 
SVM_opts.parallel_optim = true;
SVM_opts.hp_C  = 1;
SVM_opts.cv_type = "KFold"; 
SVM_opts.K = nfolds; 
SVM_opts.svm_family = 'linear'; % or 'RBF'
SVM_opts.libsvm = false;

% Define encoding boundary
encoding_W = ones(1,N_cells);

% Obtain intensities of each neuron in the population
for infoLab = info_amount
    infoLab = infoLab{1};
    delta_FR.(infoLab) = delta_FR.(infoLab)'*encoding_W;
    lambdas.(infoLab) = [mean_FR.(infoLab)-(delta_FR.(infoLab));mean_FR.(infoLab)+(delta_FR.(infoLab))]; % Hz
end

for nt=1:length(N_trials_list)
    N_trials = N_trials_list(nt);
    
    % Simulation
    R = nan(repetitions,N_cells,N_trials);
    S = nan(repetitions,N_trials);
    C = nan(repetitions,numel(angle_range),N_trials);
    
    for infoLab = info_amount
        disp('Ntrials')
        disp(N_trials)
        infoLab = infoLab{1};
        for rep = 1:repetitions
            % Simulate binary stimulus
            S(rep,:) = (binornd(1,0.5,1,N_trials));
    
            % Simulate population response
            for cell = 1:N_cells
                s0_idxs = find(S(rep,:)==0);
                [tmpR] = poissrnd(lambdas.(infoLab)(1,cell),1,numel(s0_idxs));
                R(rep,cell,s0_idxs) = tmpR;
                s1_idxs = find(S(rep,:)==1);
                [tmpR] = poissrnd(lambdas.(infoLab)(2,cell),1,numel(s1_idxs));
                R(rep,cell,s1_idxs) = tmpR;
            end
    
            for aidx = 1:numel(angle_range) % loop over encoding/decoding angles
                theta = angle_range(aidx);
                % Define 2D rotation matrix
                rotationM_2d = [cos(theta),-sin(theta); sin(theta),cos(theta)];
                % For N = 2 N-dim rotation matrix is equal to 2d
                rotationM_Nd = zeros(N_cells);
                for i = 1:N_cells/2
                    rotationM_Nd(2*i-1:2*i,2*i-1:2*i) = rotationM_2d;
                end
                decoding_W = rotationM_Nd*encoding_W';
        
                C(rep,aidx,:) = (decoding_W')*squeeze(R(rep,:,:));
                C(rep,aidx,:) = C(rep,aidx,:)>mean(squeeze(C(rep,aidx,:)));
            end
        end
        
        % Compute behavioral perf., not used (just as a check)
        for aidx = 1:numel(angle_range)
            BP.(infoLab)(aidx,:) = mean(S == squeeze(C(:,aidx,:)),2); 
        end
        
        % Decode stim and choice with SVM
        Sd = nan(repetitions,N_trials);
        Cd = nan(repetitions,numel(angle_range),N_trials);
        
        cvPartition = cvpartition(S(rep,:), 'KFold', nfolds); 
        for rep = 1:repetitions
            rep
            % Decode stim (not angle dependent)
            for cvRep = 1:nfolds
                [Sd(rep,cvPartition.test(cvRep)), ~, ~, ~, ~,W_s(rep,:),bias_s(rep)] = svm_pipeline(squeeze(R(rep,:,:))',S(rep,:)',find(cvPartition.test(cvRep)), SVM_opts); 
    
            end
            I_s.(infoLab)(rep) = cell2mat(information(Sd(rep,:)+1,S(rep,:)+1,opts,{'I'})); 
        
            for aidx = 1:numel(angle_range)
                % For each angle, decode choice
                for cvRep = 1:nfolds
                    [Cd(rep,aidx,cvPartition.test(cvRep)), ~, ~, ~, ~,W_c(rep,aidx,:),bias_c(rep,aidx)] = svm_pipeline(squeeze(R(rep,:,:))',squeeze(C(rep,aidx,:)),find(cvPartition.test(cvRep)), SVM_opts); 
                end
                I_c.(infoLab)(rep,aidx) = cell2mat(information(squeeze(Cd(rep,aidx,:))'+1,squeeze(C(rep,aidx,:))'+1,opts,{'I'})); 
                dec_SC=map_Nd_array_to_1d([Sd(rep,:);squeeze(Cd(rep,aidx,:))']);
                [ii_corr, ii_uncorr] = II(S(rep,:)+1,dec_SC,squeeze(C(rep,aidx,:))'+1,opts); 
                I_ii.(infoLab)(nt,rep,aidx,1:3) = ii_corr;
                I_ii.(infoLab)(nt,rep,aidx,4) = ii_uncorr;
                [ii_corr, ii_uncorr] = II(S(rep,:)+1,squeeze(R(rep,:,:)),squeeze(C(rep,aidx,:))'+1,opts_singlecell); 
                I_ii_singlecell.(infoLab)(nt,rep,aidx,1:3) = ii_corr;
                I_ii_singlecell.(infoLab)(nt,rep,aidx,4) = ii_uncorr;
            end
        end
    end
    
    % Average across repetitions
    slope_s = mean(-W_s(:,1)./W_s(:,2));
    slope_c = mean(-W_c(:,:,1)./W_c(:,:,2), 1);
    bias_s = mean(bias_s);
    bias_c = mean(bias_c,1);
    
end

save('ii_10-1000.mat')
%%
figure;
ax1=subplot(2,2,1);
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,1,4), 2)));hold on;
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,1,1), 2)))
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,1,2), 2)))
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,1,3), 2)))
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,1,4), 2)), squeeze(std(I_ii.high_info(:,:,1,4), [], 2))/5); hold on;
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,1,1), 2)), squeeze(std(I_ii.high_info(:,:,1,1), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,1,2), 2)), squeeze(std(I_ii.high_info(:,:,1,2), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,1,3), 2)), squeeze(std(I_ii.high_info(:,:,1,3), [], 2))/5)
legend({'Naive', "M1", "M2", "M3"})
xlabel('Number of trials per stimulus')
ylabel('II [bits]')
title('II on decoded response, angle=20deg')
% saveas(gcf, 'IIdecoded20deg101000.png')

ax2=subplot(2,2,2);
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,1,4), 2)));
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,1,1), 2)))
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,1,2), 2)))
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,1,3), 2)))
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,1,4), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,1,4), [], 2))/5);hold on;
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,1,1), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,1,1), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,1,2), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,1,2), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,1,3), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,1,3), [], 2))/5)
legend({'Naive', "M1", "M2", "M3"})
xlabel('Number of trials per stimulus')
ylabel('II [bits]')
title('II on single cell response, angle=20deg')
% saveas(gcf, 'IIsingc20deg101000.png')


ax3=subplot(2,2,3);
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,2,4), 2)));hold on;
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,2,1), 2)))
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,2,2), 2)))
% plot(N_trials_list, squeeze(mean(I_ii.high_info(:,:,2,3), 2)))
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,2,4), 2)), squeeze(std(I_ii.high_info(:,:,2,4), [], 2))/5); hold on;
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,2,1), 2)), squeeze(std(I_ii.high_info(:,:,2,1), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,2,2), 2)), squeeze(std(I_ii.high_info(:,:,2,2), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii.high_info(:,:,2,3), 2)), squeeze(std(I_ii.high_info(:,:,2,3), [], 2))/5)
legend({'Naive', "M1", "M2", "M3"})
xlabel('Number of trials per stimulus')
ylabel('II [bits]')
title('II on decoded response, angle=70deg')
% saveas(gcf, 'IIdecoded70deg101000.png')

ax4=subplot(2,2,4);
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,2,4), 2)));hold on;
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,2,1), 2)))
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,2,2), 2)))
% plot(N_trials_list, squeeze(mean(I_ii_singlecell.high_info(:,:,2,3), 2)))
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,2,4), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,2,4), [], 2))/5); hold on;
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,2,1), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,2,1), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,2,2), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,2,2), [], 2))/5)
errorbar(N_trials_list/2, squeeze(mean(I_ii_singlecell.high_info(:,:,2,3), 2)), squeeze(std(I_ii_singlecell.high_info(:,:,2,3), [], 2))/5)
legend({'Naive', "M1", "M2", "M3"})
xlabel('Number of trials per stimulus')
ylabel('II [bits]')
title('II on single cell response, angle=70deg')
linkaxes([ax1 ax2 ax3 ax4], 'y')
% saveas(gcf, 'IIsingc70deg101000.png')
saveas(gcf, ['II101000.png'])