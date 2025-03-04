classdef test_MI < matlab.unittest.TestCase
 methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()MI(), 'MI:notEnoughInput');
        end

        function testZerosInput(testCase)
            warning('off', 'all');
            S = zeros(1, 100);
            R = zeros(1, 100);
            MI_result = MI({R, S});
            assert(MI_result{1} < 0.01, 'Value should be ~zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            S = zeros(1, 100);
            S(2, 5) = NaN;
            R = zeros(1, 100);           
            opts.NaN_handling = 'error';
            opts.suppressWarnings = true;
            testCase.verifyError(@()MI({R, S}, opts), 'checkInputs:NaNDetected'); 
        end

        function testRandiValues(testCase)
            warning('off', 'all');
            rng("default")
            S = randi([1, 3], 1, 100);
            R = randi([1, 3], 1, 100);           
            MI_result = MI({R, S});
            assert(MI_result{1} < 0.02, 'Value should be ~zero.');
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
            R = [];
            S = [];
            C = [];
            for i = 1:nStimuli
                for j = 1:nTrials
                    R = [R sum(poisson_spike_gen(t_trial, stimuli(i,:), 0))];
                    S = [S i];
                end
            end
            opts.bias = 'plugin';
            opts.bin_method = {'eqpop'};
            opts.suppressWarnings  = true;
            opts.n_bins = {2, 3};
            MI_out = MI({R, S},opts);
            assert((MI_out{1} - 0.6205) < 0.01);
            opts.bias = 'invalid';
            testCase.verifyError(@()MI({R, S},opts), 'Correction:UndefinedFunction');
            opts.bias = 'plugin';
            opts.bin_method  = {'invalid'};
            testCase.verifyError(@()MI({R, S},opts), 'Binning:UndefinedFunction');
        end
    end

end