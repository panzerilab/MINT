classdef test_cMI < matlab.unittest.TestCase
 methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()cMI(), 'cMI:notEnoughInput');
        end

        function testZerosInput(testCase)
            warning('off', 'all');
            S = zeros(1, 100);
            R1 = zeros(1, 100);
            R2 = zeros(1, 100);
            cMI_result = cMI({R1, S, R2});
            assert(cMI_result{1} < 0.01, 'Value should be ~zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            S = zeros(1, 100);
            S(2, 5) = NaN;
            R1 = zeros(1, 100);
            R2 = zeros(1, 100);       
            opts.NaN_handling = 'error';
            opts.supressWarnings = true;
            testCase.verifyError(@()cMI({R1, S, R2}, opts), 'checkInputs:NaNDetected'); 
        end

        function testRandiValues(testCase)
            warning('off', 'all');
            rng("default")
            S = randi([1, 2], 1, 100);
            R1 = randi([1, 2], 1, 100);
            R2 = randi([1, 2], 1, 100);
            cMI_result = cMI({R1, S, R2});
            assert(cMI_result{1} < 0.075, 'Value should be ~zero.');
            warning('on', 'all');
        end
       
        function testValid_InvalidInputs(testCase)
            rng("default")
            dt = 1/50;
            trialendtime = 1;
            t_trial = 0:dt:trialendtime;
            nStimuli = 2;
            nTrials = 20;
            mu = .2;  
            sigmaT = [.05,.05];
            rate = [10 50];  
            stimuli = nan(nStimuli,length(t_trial));
            for i = 1:nStimuli
                signal = normpdf(t_trial,mu,sigmaT(i));
                stimuli(i,:) =  rate(i) * signal / max(signal); 
            end
            R1 = [];
            S = [];
            for i = 1:nStimuli
                for j = 1:nTrials
                    R1 = [R1 sum(poisson_spike_gen(t_trial, stimuli(i,:), 0))];
                    S = [S i];
                end
            end
            R2 = randi([1,3],1, 40);
            opts.bias = 'naive';
            opts.bin_method = {'eqpop', 'none', 'eqpop'};
            opts.supressWarnings  = true;
            opts.n_bins = {2, 3};
            cMI_out = cMI({R1, S, R2},opts);
            assert((cMI_out{1} - 0.6833) < 0.01);
            opts.bias = 'invalid';
            testCase.verifyError(@()cMI({R1, S, R2},opts), 'Correction:UndefinedFunction');
            opts.bias = 'naive';
            opts.bin_method  = {'invalid'};
            testCase.verifyError(@()cMI({R1, S, R2},opts), 'Binning:UndefinedFunction');
        end
    end

end