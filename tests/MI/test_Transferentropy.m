classdef test_Transferentropy < matlab.unittest.TestCase

    methods (Test)
        function testEmptyInput(testCase)
            testCase.verifyError(@()transferentropy(), 'transferentropy:notEnoughInput');
        end

        function testTransferentropyRandomInput(testCase)
            rng('default');
            X = rand(1,100, 50);
            Y = rand(1,100, 50);
            opts.taux = -1;
            opts.tauy = -1;
            opts.bias = 'naive';
            opts.method = 'dr';
            opts.bin_methodX = 'eqpop';
            opts.bin_methodY = 'eqpop';
            opts.n_binsX = 3;
            opts.n_binsY = 3;
            outputList = {'TE'};
            opts.verbose = false;
            opts.supressWarnings = true;
            TE = transferentropy(X, Y, opts, outputList);
            tol = 0.01;
            TE{1, 1} <= tol;
        end

        function testTransferentropyValidInput(testCase)
            rng('default');
            nTrials = 100;
            nTimepoints = 20;
            noiseStr = 0.1;
            S = binornd(1,0.5,1,nTrials);
            S_reshaped = reshape(S, [1, 1, nTrials]);
            S_extended = repmat(S_reshaped, [1, nTimepoints, 1]);
            X =  S_extended + noiseStr*randn(1,nTimepoints, nTrials);
            Y = X + noiseStr * randn(1, nTimepoints, nTrials);
            opts.method = 'dr';                                             
            opts.bias   = 'naive'; 
            opts.bin_methodX = 'none';
            opts.bin_methodY = 'none';
            opts.singleTimepoint = true;
            opts.supressWarnings = true;
            opts.verbose = 0;
            opts.taux = -5;
            opts.tauy = -5;
            outputList = {'TE'};
            opts.verbose = false;
            TE = transferentropy(X, Y, opts, outputList);
            tol = 0.01;
            abs(TE{1} - 0.3759) <= tol;
        end


        function testTransferentropyNaNInput(testCase)
            rng('default');
            X = rand(100, 3);
            Y = rand(100, 3);
            Y(2,2) = NaN;
            opts.taux = -1;
            opts.supressWarnings = true;
            opts.tauy = -1;
            opts.bias = 'naive';
            opts.method = 'dr';
            outputList = {'TE'};
            opts.verbose = false;
            testCase.verifyError(@() transferentropy(X, Y, opts, outputList), 'sanitycheck:NaNInput');
        end

        function testTransferentropyInvalidInput(testCase)
            rng('default');
            X = rand(100, 3);
            Y = rand(100, 3);
            opts.supressWarnings = true;
            opts.taux = -1;
            opts.tauy = -1;
            opts.bias = 'invalid';
            opts.method = 'dr';
            outputList = {'TE'};
            opts.verbose = false;
            testCase.verifyError(@() transferentropy(X, Y, opts, outputList), 'transferentropy:InvalidBiasMethod');
        end

        function testTransferentropyMissingInput(testCase)
            rng('default');
            X = rand(100, 3);
            Y = rand(100, 3);
            opts.supressWarnings = true;
            opts.taux = -1;
            opts.tauy = -1;
            opts.bias = 'naive';
            outputList = {'TE'};
            opts.verbose = false;
            testCase.verifyError(@() transferentropy(X, Y, opts, outputList), 'transferentropy:MissingInput');
        end
    end
end