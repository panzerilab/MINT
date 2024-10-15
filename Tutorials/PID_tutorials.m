% -------------------------------------------------------------------------
%                       TUTORIAL PID
% -------------------------------------------------------------------------

%% Partial Information Decomposition (PID)

dt = 1/50;
trialendtime = 1;
t_trial = 0:dt:trialendtime;
nStimuli = 2;
nTrials = 100;

% GENERATE STIMULUS
mu = .2;                                            % mu for the gaussian
sigmaT = [.05,.05];                                 % sigma for each stimulus
rate = [10 50];                                     % peak rate for each stimulus
stimuli = nan(nStimuli,length(t_trial));
for i = 1:nStimuli
    signal = normpdf(t_trial,mu,sigmaT(i));
    stimuli(i,:) =  rate(i) * signal / max(signal); %normalize for peak
end

% GENERATE RESPONSE AND CHOICE
R = [];
S = [];
C = [];
for i = 1:nStimuli
    for j = 1:nTrials
        R = [R sum(poisson_spike_gen(t_trial, stimuli(i,:), 0))];
        S = [S i];
    end
end
C = zeros(size(R));
C(R > 2) = 1;


% OPTIONS TO CALCULATE PID
opts_PID.bias = 'naive';                               % bias correction    
opts_PID.bin_methodX = 'eqspace';                  % binning
opts_PID.n_binsX = 2;
opts_PID.bin_methodY = 'eqpop';
opts_PID.n_binsY = 3;
opts_PID.shuff = 10;                                % bootstrapping
opts_PID.shuff_variables = {'Y'};
opts_PID.shuff_type = {'all'};
opts_PID.shuff_bias = {'naive'};
opts_PID.redundancy_measure = "I_BROJA";              % redundancy mesure method
opts_PID.xtrp = 4;

% CALCULATE PID
[PID_v, uncorrected] = PID({R; S; R}, C, opts_PID);

fprintf("PID_SR_C shared =  %.3f\n", PID_v(1))
fprintf("PID_SR_C uniqueR =  %.3f\n", PID_v(2))
fprintf("PID_SR_C uniqueS =  %.3f\n", PID_v(3))
fprintf("PID_SR_C complementary  =  %.3f\n", PID_v(4))

%% FEATURE-SPECIFIC INFORMATION TRANSFER (FIT)

% GENERATE STIMULUS
nTrials = 10000;
noiseStr = 0.5;
S = binornd(1,0.5,1,nTrials); % Binomial stimulus

% GENERATE NEURAL RESPONSE  S -> X -> Z -> Y
X_past = S + noiseStr*randn(1,nTrials); 
Z_past = X_past + noiseStr*randn(1,nTrials);
Y_pres = Z_past + noiseStr*randn(1,nTrials);
Y_past = noiseStr*randn(1,nTrials);

% OPTIONS TO CALCULATE FIT
opts_FIT.verbose = false;
opts_FIT.method = "dr";
opts_FIT.bias = 'qe';
opts_FIT.bin_methodX = 'eqpop';
opts_FIT.bin_methodY = 'eqpop';
opts_FIT.n_binsX = 5;
opts_FIT.n_binsY = 4;
opts_FIT.bin_methodS = 'none';
opts_FIT.shuff = 0;
opts_FIT.shuff_variables = {'hX'};
opts_FIT.shuff_type = {'Sconditioned'};
opts_FIT.shuff_bias = {'naive'};
opts_FIT.parallel = 1;

% COMPUTE FIT (S:X->Y)
% [FIT_v, shuff_v, FIT_v_uncorrected] =FIT(S, X_past, Y_past, Y_pres, opts_FIT);
disp(['Measured FIT = ', num2str(FIT_v.biased.value) , ' pval = ', num2str(mean(FIT_v.biased.shuff{1} >= FIT_v.biased.value))]);