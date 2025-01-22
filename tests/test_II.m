classdef test_II < matlab.unittest.TestCase

    methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()II(), 'II:notEnoughInput');
        end

        function testZerosInput(testCase)
            warning('off', 'all');
            S = zeros(1, 100);
            R = zeros(1, 100);
            C = zeros(1, 100);
            II_result = II({S, R, C});
            assert(II_result{1} < 0.01, 'Value should be ~zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            S = rand(3, 100);
            R = rand(3, 100);
            C = rand(3, 100);
            S(2, 5) = NaN;
            opts.NaN_handling = 'error';
            opts.suppressWarnings = true;
            testCase.verifyError(@()II({S, R, C}, opts), 'checkInputs:NaNDetected'); 
        end

        function testRandiValues(testCase)
            warning('off', 'all');
            rng("default")
            S = randi([1, 3], 1, 100);
            R = randi([1, 3], 1, 100);
            C = randi([1, 3], 1, 100);
            II_result = II({S, R, C});
            assert(II_result{1} < 0.02, 'Value should be ~zero.');
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
            C = zeros(size(R));
            C(R > 2) = 1;
            opts.bias = 'plugin';
            opts.bin_method = {'eqpop'};
            opts.suppressWarnings  = true;
            opts.n_bins = {2, 3, 2};
            II_out = II({S, R, C},opts);
            assert((II_out{1} - 0.5969) < 0.01);
            opts.bias = 'invalid';
            testCase.verifyError(@()II({S, R, C},opts), 'Correction:UndefinedFunction');
            opts.bias = 'plugin';
            opts.bin_method  = {'invalid'};
            testCase.verifyError(@()II({S, R, C},opts), 'Binning:UndefinedFunction');
        end
    end

end