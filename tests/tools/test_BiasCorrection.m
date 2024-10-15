classdef test_BiasCorrection < matlab.unittest.TestCase

    methods (Test) 
        function bias_convergence(testCase)
            rng('default');
            dt = 1/500;
            trialendtime = 0.4;
            t_trial = 0:dt:trialendtime;
            nStimuli = 2;

            shuffopts.bias = 'shuffSub';
            shuffopts.bin_method = {'eqpop', 'none'};
            shuffopts.n_bins = {3, 0};
            shuffopts.shuff = 5;
            shuffopts.compute_nulldist = false;
            shuffopts.parallel = 0;

            qeopts = shuffopts;
            qeopts.bias = 'qe';
            qeopts.xtrp = 5;
            qeopts.shuff = 0;

            shopts = shuffopts;
            shopts.bias = 'shuffCorr';
            shopts.xtrp = 0;
            shopts.shuff = 0;
            nTrials_list = [64, 256, 512, 1000];
            % generate response to stimuli
            maxrate = [50 20];
            rate = nan(nStimuli,length(t_trial));
            for i = 1:nStimuli
                signal = normpdf(t_trial,0.2,0.05);
                rate(i,:) =  maxrate(i) * signal / max(signal);
            end

            MI_naive = zeros(1,length(nTrials_list));
            MI_shuff = zeros(1,length(nTrials_list));
            MI_qe    = zeros(1,length(nTrials_list));
            MI_sh    = zeros(1,length(nTrials_list));

            parfor tidx = 1:length(nTrials_list)
                % generate responses (spike count)
                R = [];
                S = [];
                nTrials = nTrials_list(tidx);
                for i = 1:nStimuli
                    for j = 1:nTrials
                        spike_train = poisson_spike_gen(t_trial, rate(i,:), 0);
                        R = [R sum(spike_train)];
                        S = [S i];
                    end
                end
                
                [mi_shuff, mi_naive] = MI({R;S}, {'I(A;B)'},shuffopts);
                mi_qe                = MI({R;S}, {'I(A;B)'},qeopts);
                mi_sh                = MI({R;S}, {'I(A;B)'},shopts);

                MI_naive(tidx) = mi_naive;
                MI_shuff(tidx) = mi_shuff;
                MI_qe(tidx)    = mi_qe;
                MI_sh(tidx)    = mi_sh;

            end

            figure;
            semilogx(nTrials_list, MI_naive)
            semilogx(nTrials_list, MI_shuff)
            semilogx(nTrials_list, MI_qe   )
            semilogx(nTrials_list, MI_sh   )

            legend({'MI_naive', 'MI_shuff', 'MI_qe', 'MI_sh'})

            % testCase.verifyError(@() information(R, S, opts, {'I'}), 'sanitycheck:NaNInput');        
            % functionCall = @() gaussian_bias();
            % verifyError(testCase, functionCall, 'MATLAB:minrhs', 'Expected MATLAB:minrhs error.');
        end
        
    end
end