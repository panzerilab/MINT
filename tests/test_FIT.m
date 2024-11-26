classdef test_FIT < matlab.unittest.TestCase

    methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()FIT(), 'FIT:notEnoughInput');
        end

        function testZerosInput(testCase)
            warning('off', 'all');
            S = zeros(1, 100);
            X = zeros(1, 100, 100);           
            Y = zeros(1, 100, 100);
            FIT_result = FIT({X, Y, S});
            assert(FIT_result{1} == 0, 'Value should be zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            S = rand(3, 100);
            X = rand(3, 100, 100);
            Y = rand(3, 100, 100);           
            Y(2, 5) = NaN;
            opts.NaN_handling = 'error';
            opts.suppressWarnings = true;
            testCase.verifyError(@()FIT({X, Y, S}, opts), 'checkInputs:NaNDetected'); 
        end

        function testValidInput(testCase)
            rng("default")
            nTrials = 100;
            nTimepoints = 20;
            noiseStr = 0.1;
            S = binornd(1,0.5,1,nTrials);
            S_reshaped = reshape(S, [1, 1, nTrials]);
            S_extended = repmat(S_reshaped, [1, nTimepoints, 1]);
            X =  S_extended + noiseStr*randn(1,nTimepoints, nTrials);
            Y = noiseStr * randn(1, nTimepoints, nTrials);
            Y(:,1:(nTimepoints/2),:) = noiseStr*randn(1,nTimepoints/2, nTrials);
            Y(:, (nTimepoints/2 + 1):nTimepoints, :) = X(:, (nTimepoints/2 + 1):nTimepoints, :) + noiseStr * randn(1, nTimepoints/2, nTrials);
            opts.suppressWarnings = true;            
            opts.bin_method = {'eqpop','eqpop', 'none'};
            opts.n_bins = {3, 3};
            opts.tau = {12};
            opts.tpres = {20};
            FIT_result = FIT({X, Y, S}, opts);
            assert((FIT_result{1} - 0.4260) < 0.01);
        end

        function testValidInput_qe(testCase)
             rng("default")
            nTrials = 100;
            nTimepoints = 20;
            noiseStr = 0.1;
            S = binornd(1,0.5,1,nTrials);
            S_reshaped = reshape(S, [1, 1, nTrials]);
            S_extended = repmat(S_reshaped, [1, nTimepoints, 1]);
            X =  S_extended + noiseStr*randn(1,nTimepoints, nTrials);
            Y = noiseStr * randn(1, nTimepoints, nTrials);
            Y(:,1:(nTimepoints/2),:) = noiseStr*randn(1,nTimepoints/2, nTrials);
            Y(:, (nTimepoints/2 + 1):nTimepoints, :) = X(:, (nTimepoints/2 + 1):nTimepoints, :) + noiseStr * randn(1, nTimepoints/2, nTrials);          
            opts.suppressWarnings = true;
            opts.bias = 'qe_shuffSub';
            opts.xtrp = 10;
            opts.bin_method = {'eqpop','eqpop', 'none'};
            opts.n_bins = {3, 3};
            opts.tau = {12};
            opts.tpres = {20};
            FIT_result = FIT({X, Y, S}, opts);
            assert((FIT_result{1} - 0.4367) < 0.01);
        end


        function testInvalidInput(testCase)
           rng('default');
            warning('off', 'all');
            S = zeros(1, 100);
            X = zeros(1, 100, 100);           
            Y = zeros(1, 100, 100);
            opts.bias = 'invalid';
            opts.xtrp= 20;
            opts.suppressWarnings = true;
            opts.bin_method = {'eqspace', 'eqpop', 'none'};
            opts.n_bins = {2, 2};
            testCase.verifyError(@()FIT({X, Y, S}, opts), 'Correction:UndefinedFunction');

            opts.bias = 'naive';
            opts.bin_method  = {'invalid'};
            testCase.verifyError(@()FIT({X, Y, S}, opts), 'Binning:UndefinedFunction');
            warning('on', 'all');
        end
    end
end 