%% TUTORIAL MINT
clc, clear
close all;

% -------------------------------------------------------------------------
%                            Simulate Neural Response
% -------------------------------------------------------------------------

% Parameters 
nTrials = 800;
t = 0:0.5:10000; %ms
dt = t(2)-t(1);
nTimeSteps = length(t);

% Populations
nNeurons_1 = 8;
nNeurons_2 = 8;
nNeurons = nNeurons_1+nNeurons_2;
selectivity_1 = [zeros(1, nNeurons_1/2), ones(1, nNeurons_1/2)+2];
selectivity_2 = [zeros(1, nNeurons_2/2), ones(1, nNeurons_2/2)+2];
selectivity = [selectivity_1, selectivity_2];
isExcitatory = (mod(1:nNeurons, 2) == 1)';
scaleFactor = 1 / sqrt(nNeurons);  
excitatoryWeight = 8*scaleFactor;
inhibitoryWeight = -4*scaleFactor;
p = 1;
connections = zeros(nNeurons, nNeurons);
connections(1:nNeurons_1, 1:nNeurons_1) = rand(nNeurons_1, nNeurons_1) < p;
connections(1:nNeurons_1, nNeurons_1+1:nNeurons) = rand(nNeurons_1, nNeurons_2) < p;
connections(nNeurons_1+1:nNeurons, 1:nNeurons_1) = rand(nNeurons_2, nNeurons_1) < p;
connections(nNeurons_1+1:nNeurons, nNeurons_1+1:nNeurons) = rand(nNeurons_2, nNeurons_2) < p;

baseCorrelation = 0.1; 
correlationMatrix = baseCorrelation * ones(nNeurons);
maxCorrelation = 0.1; 
minCorrelation = -0.1;

for i = 1:nNeurons
    for j = 1:nNeurons
        if i == j
            correlationMatrix(i, j) = 1; 
        else
            diff = abs(selectivity(i) - selectivity(j));
            if diff == 0
                correlationMatrix(i, j) = maxCorrelation; 
            else
                correlationMatrix(i, j) = minCorrelation; 
            end
        end
    end
end

epsilon = 1e-4;
correlationMatrix = correlationMatrix + epsilon * eye(nNeurons);
L = chol(correlationMatrix, 'lower');

w = zeros(nNeurons, nNeurons);
% Weights Population 1
for i = 1:nNeurons_1
    for j = 1:nNeurons_1
        if connections(i, j)
            diff = abs(selectivity(i) - selectivity(j));
            if isExcitatory(j)
                if diff == 0
                    w(i, j) = excitatoryWeight;
                end
            else
                if diff ~= 0 && i ~= j && isExcitatory(i)
                    w(i, j) = inhibitoryWeight;
                end
            end
        end
    end
end

% Weights Population 1 -> Population 2
for i = nNeurons_1+1:nNeurons
     for j = 1:nNeurons_1
        if connections(i, j)
            diff = abs(selectivity(i) - selectivity(j));
            if isExcitatory(j)
                if diff == 0
                    w(i, j) = excitatoryWeight*10;
                end
            else
                if diff ~= 0 && isExcitatory(i)
                    w(i, j) = inhibitoryWeight*2;
                end
            end
        end
     end
end

% Weights Population 2
for i = nNeurons_1+1:nNeurons
    for j = nNeurons_1+1:nNeurons
        if connections(i, j)
            diff = abs(selectivity(i) - selectivity(j));
            if isExcitatory(j)
                if diff == 0
                    w(i, j) = excitatoryWeight;
                end
            else
                if diff ~= 0 && i ~= j && isExcitatory(i)
                    w(i, j) = inhibitoryWeight;
                end
            end
        end
    end
end

% figure;
% imagesc(w);
% colorbar;
% axis square;
% xticks(1:nNeurons);
% yticks(1:nNeurons);
% line([0.5, nNeurons+0.5], [nNeurons_1+0.5, nNeurons_1+0.5], 'Color', 'black', 'LineWidth', 2);
% line([nNeurons_1+0.5, nNeurons_1+0.5], [0.5, nNeurons+0.5], 'Color', 'black', 'LineWidth', 2);


% Generate Stimulus
stimulusDirection = [zeros(1,nTrials/4), ones(1,nTrials/4), ones(1,nTrials/4)+1, ones(1,nTrials/4)+2];
inputIntensity_direction = 5;
color = [zeros(1,nTrials/8), ones(1,nTrials/8)];
stimulusColor = [color, color, color, color];

% Network Parameter
V_0 = -70;
V_th = -50;
tau = 30;
V = zeros(nNeurons, length(t),nTrials);
M = zeros(nNeurons, length(t),nTrials);
c = 30;

for trial = 1:nTrials

    if stimulusColor(trial) == 0
        inputColor = 2;
    else 
        inputColor = 1;
    end 

    inputs = zeros(nNeurons, nTimeSteps);
    for i = 1:nNeurons
        if i <= nNeurons/2
            diff = selectivity(i) - stimulusDirection(trial);
            if diff == 0
                inputs(i, :) = (inputIntensity_direction*inputColor * (sin(2 * pi * 5 * t / max(t)-pi/2) + 1) / 2) + 3 * randn(1, length(t));
            elseif abs(diff) == 1
                inputs(i, :) = (inputIntensity_direction*inputColor/1.5 * (sin(2 * pi * 5 * t / max(t)-pi/2) + 1) / 2) + 3 * randn(1, length(t));
            else
                inputs(i, :) = normrnd(0, 15, 1, length(t))/sqrt(dt)*inputColor;
            end
        else
            inputs(i, :) = normrnd(0, 10, 1, length(t))/sqrt(dt);
        end
    end

    % noise
    uncorrelatedNoise = normrnd(0, 75, nNeurons, nTimeSteps);
    correlatedNoise = L * uncorrelatedNoise;
 
    V(:,1, trial)=V_0;
    for i = 2:length(t)
        dV = -(V(:,i-1, trial) - V_0) + c*w*M(:,i-1,trial)+ inputs(:, i) + correlatedNoise(:,i);
        dV = dV/tau*dt;
        V(:,i, trial) = V(:,i-1, trial)+dV;
        for n = 1:nNeurons
            if V(n,i, trial) > V_th
                V(n,i-1, trial) = 0;
                M(n,i, trial) = 1;
                V(n,i, trial) = V_0;
            end
        end
    end
end

V_c = V;
V(V >= -0.01) = 1;
V(V < 0) = 0;

figure
for n = 1:nNeurons
    plot(t, n*V(n, :, 2),'.', 'Color', 'red');  
    hold on
    plot(t, n*V(n, :, (nTrials/2)+2),'.', 'Color', 'blue');  
end

%%
R = squeeze(sum(V,2));
choice = zeros(1, nTrials);
for t = 1:nTrials
    if sum(R(9:12, t)) > sum(R(13:16, t))
        choice(t) = 1;
    else 
        choice(t) = 2;
    end 
end 

R_mean_0 = mean(R(:, 1:nTrials/4),2);
R_mean_1 = mean(R(:, (nTrials/4)+1:nTrials/2), 2);
R_mean_2 = mean(R(:, (nTrials/2)+1:(3*nTrials)/4), 2);
R_mean_3 = mean(R(:, ((3*nTrials)/4)+1:nTrials), 2);

% -------------------------------------------------------------------------
%                            MUTUAL INFORMATION
% -------------------------------------------------------------------------
%%
opts_MI.bias = 'naive';
opts_MI.method = 'dr';
opts_MI.outputsList = {'HX', 'HXY'};
opts_MI.bin_methodX = 'eqspace';
opts_MI.n_binsX = 3;
opts_MI.bin_methodY = 'none';
opts_MI.verbose = 0;

Info = zeros(1, nNeurons);
for n = 1:nNeurons
    outputs = information(R(n,:), [stimulusDirection; stimulusColor], opts_MI, {'I'});
    Info(n) = outputs{1};
end

figure 
plot(Info);

% -------------------------------------------------------------------------
%                            INFORMATION BREAKDOWN
% -------------------------------------------------------------------------
%%
info_bdw_terms = {'Joint','ILIN','ISS','ICI','ICD'};
nPairs = nNeurons*(nNeurons-1)/2;
MI_breakdown = struct(); 

for bdwIdx = 1:numel(info_bdw_terms)
    bdwLab = info_bdw_terms{bdwIdx};
    MI_breakdown.(bdwLab) = zeros(nNeurons, nNeurons);
end
    

pairIdx = 0;
for neuron1 = 1:nNeurons
    for neuron2 = neuron1+1:nNeurons
        jointResp = [R(neuron1,:); R(neuron2,:)];
        infoBdw = information(jointResp, stimulusDirection, opts_MI, {'I','ILIN','ISS','ICI','ICD'});
        for bdwIdx = 1:numel(info_bdw_terms)
            bdwLab = info_bdw_terms{bdwIdx};           
            MI_breakdown.(bdwLab)(neuron1, neuron2) = infoBdw{bdwIdx};
        end
    end
end

% -------------------------------------------------------------------------
%                 Population Information with SVM Decoder
% -------------------------------------------------------------------------

% Set SVM Parameters
opts_svm.optimize_params = false;   
opts_svm.cv_type = "KFold"; 
opts_svm.K = 5; 
opts_svm.libsvm = false; 

cvPartition = cvpartition(stimulusDirection, opts_svm.cv_type, opts_svm.K);

% Linear SVM
info_out_all = [];
for i = 1:cvPartition.NumTestSets
    test_idxs = find(cvPartition.test(i));
    opts_svm.svm_family = 'linear';
    [PredLabels] = svm_pipeline(R',stimulusDirection',test_idxs, opts_svm);
    info_out = information(PredLabels', stimulusDirection(test_idxs), opts_MI, {'I'});
    info_out_all(end+1) = info_out{1};
end
info_out_mean = mean(info_out_all);
MI.linear.mean = info_out_mean;
MI.linear.all = info_out_all;

% RBF SVM
info_out_all = [];
for i = 1:cvPartition.NumTestSets
    test_idxs = find(cvPartition.test(i));
    opts_svm.svm_family = 'RBF';
    [PredLabels] = svm_pipeline(R',stimulusDirection',test_idxs, opts_svm);
    info_out = information(PredLabels', stimulusDirection(test_idxs), opts_MI, {'I'});
    info_out_all(end+1) = info_out{1};
end
info_out_mean = mean(info_out_all);
MI.RBF.mean = info_out_mean;
MI.RBF.all = info_out_all;


% -------------------------------------------------------------------------
%                          INTERSECTION INFORMATION
% -------------------------------------------------------------------------

opts_II.bias = 'naive';
opts_II.method = 'dr';
opts_II.bin_methodS = 'none';
opts_II.bin_methodC = 'none';
opts_II.bin_methodR = 'none';
opts_II.verbose = 0;

II_inf = zeros(1, nNeurons);

for i = 1:nNeurons
        II_inf(i) = II(stimulusDirection, R(i,:), choice, opts_II);      
end

%%

% -------------------------------------------------------------------------
%                            TRANSFER ENTROPY
% -------------------------------------------------------------------------
opts_TE.bias = 'naive';
opts_TE.method = 'dr';
outp_TEutsList = {'HX', 'HXY'};
opts_TE.bin_methodX = 'eqspace';
opts_TE.n_binsX = 2;
opts_TE.bin_methodY = 'eqspace';
opts_TE.n_binsY = 2;
opts_TE.verbose = 0;
opts_TE.taux   = [-1];
opts_TE.tauy   = [-1];

TE = zeros(nNeurons, nNeurons);

for i = 1:nNeurons
    for j = 1:nNeurons
        outputs_1 = transferentropy(squeeze(V(i,:,:)), squeeze(V(j,:,:)), opts_TE, {'TE'});
        TE(i,j) = outputs_1{1};
    end
end

figure;
hold on;
xlabel('Connection');
ylabel('Transferentropie');

connected_TE = [];
unconnected_TE = [];

for i = 1:nNeurons
    for j = 1:nNeurons
        if i ~= j
            if w(i,j) > 0
                connected_TE(end+1) = TE(j, i);
                plot(1, TE(j,i), 'bo');
            else
                unconnected_TE(end+1) = TE(j, i);
                plot(0, TE(j,i), 'rx');
            end
        end
    end
end

mean_connected_TE = mean(connected_TE);
mean_unconnected_TE = mean(unconnected_TE);
plot([0 1], [mean_unconnected_TE mean_connected_TE], 'g-', 'LineWidth', 1);

R_coeff = corrcoef([zeros(size(unconnected_TE)), ones(size(connected_TE))], [unconnected_TE, connected_TE]);
disp(['R: ', num2str(R_coeff(1,2))]);

hold off;





%%

% -------------------------------------------------------------------------
%                                 PID
% -------------------------------------------------------------------------
opts_PID.bias = 'naive';
opts_PID.method = 'dr';
outp_PIDutsList = {'HX', 'HXY'};
opts_PID.bin_methodX = 'eqspace';
opts_PID.n_binsX = 2;
opts_PID.bin_methodY = 'eqspace';
opts_PID.n_binsY = 2;
opts_PID.verbose = 0;

PID_value = PID(stimulusDirection, R(1,:), R(5,:), opts);

