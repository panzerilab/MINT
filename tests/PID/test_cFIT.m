%% 
classdef test_cFIT < matlab.unittest.TestCase

    methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()cFIT(), 'cFIT:notEnoughInput');
        end

        function testZerosInput(testCase)
            warning('off', 'all');
            S = zeros(3, 100);
            X = zeros(3, 100, 100);           
            Y = zeros(3, 100, 100);
            Z = zeros(3, 100, 100);               
            cFIT_result = cFIT(S, X, Y, Z);
            assert(cFIT_result == 0, 'Value should be zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            S = rand(3, 100);
            X = rand(3, 100, 100);
            Y = rand(3, 100, 100);
            Z = rand(3, 100, 100);
            Y(2, 5) = NaN;
            testCase.verifyError(@()cFIT(S, X, Y, Z), 'sanitycheck:NaNInput');
        end

        function testRandiValues(testCase)
            warning('off', 'all');
            rng("default")
            S  = randi([1, 2], 1, 100);
            X = randi([1, 2], 1, 100, 100);
            Y = randi([1, 2], 1, 100, 100);         
            Z = randi([1, 2], 1, 100, 100);
            cFIT_result = cFIT(S, X, Y, Z);
            assert(cFIT_result == 0);
            warning('on', 'all');
        end

        function testRandValues(testCase)
            rng("default")
            S  = rand(1, 100);
            X = rand(1, 100, 100);
            Y = rand(1, 100, 100);
            Z = rand(1, 100, 100);
            cFIT_result = cFIT(S, X, Y, Z);
            assert(cFIT_result < 0.01);
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
            Z = noiseStr * randn(1, nTimepoints, nTrials);
            Z(:,1:(nTimepoints/2),:) = noiseStr*randn(1,nTimepoints/2, nTrials);
            Z(:, (nTimepoints/2 + 1):nTimepoints, :) = X(:, (nTimepoints/2 + 1):nTimepoints, :) + noiseStr * randn(1, nTimepoints/2, nTrials);
            Y = noiseStr * randn(1, nTimepoints, nTrials);
            Y(:,1:(nTimepoints/2),:) = noiseStr*randn(1,nTimepoints/2, nTrials);
            Y(:, (nTimepoints/2 + 1):nTimepoints, :) = Z(:, (nTimepoints/2 + 1):nTimepoints, :) + noiseStr * randn(1, nTimepoints/2, nTrials);
            opts.verbose = false;
            opts.method = "dr";
            opts.bias = 'naive';
            opts.bin_methodX = 'eqpop';
            opts.bin_methodY = 'eqpop';
            opts.bin_methodZ = 'eqpop';
            opts.n_binsX = 5;
            opts.n_binsY = 5;
            opts.n_binsZ = 5;
            opts.taux = -12;
            opts.tauz = -12;
            opts.tauy = -12;
            opts.tpres = 20;
            opts.bin_methodS = 'none';
            cFIT_result = cFIT(S, X, Y, Z, opts);
            assert((cFIT_result - 0.5169) < 0.01);
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
            Z = noiseStr * randn(1, nTimepoints, nTrials);
            Z(:,1:(nTimepoints/2),:) = noiseStr*randn(1,nTimepoints/2, nTrials);
            Z(:, (nTimepoints/2 + 1):nTimepoints, :) = X(:, (nTimepoints/2 + 1):nTimepoints, :) + noiseStr * randn(1, nTimepoints/2, nTrials);
            Y = noiseStr * randn(1, nTimepoints, nTrials);
            Y(:,1:(nTimepoints/2),:) = noiseStr*randn(1,nTimepoints/2, nTrials);
            Y(:, (nTimepoints/2 + 1):nTimepoints, :) = Z(:, (nTimepoints/2 + 1):nTimepoints, :) + noiseStr * randn(1, nTimepoints/2, nTrials);
            opts.verbose = false;
            opts.method = "dr";
            opts.bias = 'qe';
            opts.bin_methodX = 'eqpop';
            opts.bin_methodY = 'eqpop';
            opts.bin_methodZ = 'eqpop';
            opts.n_binsX = 5;
            opts.n_binsY = 5;
            opts.n_binsZ = 5;
            opts.taux = -12;
            opts.tauz = -12;
            opts.tauy = -12;
            opts.tpres = 20;
            opts.bin_methodS = 'none';
            cFIT_result = cFIT(S, X, Y, Z, opts);
            assert((cFIT_result - 0.6541) < 0.01);
        end

        function testInvalidBias(testCase)
            rng("default")
            S  = randi([1, 2], 1, 100);
            X = randi([1, 2], 1, 100, 100);
            Y = randi([1, 2], 1, 100, 100);         
            Z = randi([1, 2], 1, 100, 100);
            opts.verbose = false;
            opts.bias = 'invalid';
            opts.bin_methodX = 'eqpop';
            opts.bin_methodY = 'eqpop';
            opts.n_binsX = 5;
            opts.n_binsY = 5;
            opts.n_binsZ = 3;
            opts.bin_methodS = 'none';
            testCase.verifyError(@()cFIT(S, X, Y, Z, opts), 'cFIT:UndefinedBiasCorrectionMethod');
        end

        function testInvalidBinning(testCase)
            rng("default")
            S  = randi([1, 2], 1, 100);
            X = randi([1, 2], 1, 100, 100);
            Y = randi([1, 2], 1, 100, 100);         
            Z = randi([1, 2], 1, 100, 100);
            opts.verbose = false;
            opts.bias = 'naive';
            opts.bin_methodX = 'invalid';
            opts.bin_methodY = 'eqpop';
            opts.n_binsX = 5;
            opts.n_binsY = 5;
            opts.n_binsZ = 3;
            opts.bin_methodS = 'none';
            testCase.verifyError(@()cFIT(S, X, Y, Z, opts), 'cFIT:InvalidBinningMethod');
        end
    end

end