classdef test_Entropy < matlab.unittest.TestCase

    methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()entropy(), 'entropy:notEnoughInput');
        end
        % 
        % function testValidInputs(testCase)
        %     rng('default');
        %     nTrials = 400;
        %     nDimensionsX = 3;
        %     nDimensionsY = 2;
        %     X = randn(nDimensionsX, nTrials);
        %     Y = randi([1, 5], nDimensionsY, nTrials);
        %     opts.n_binsX = 5;
        %     opts.n_binsY = 5;
        %     opts.verbose = false;
        %     opts.bias = 'gsb';
        %     opts.method = 'gs';
        %     outputsList = {'HX', 'HXY'};
        %     entropies = entropy(X, Y, opts, outputsList);
        %     expected_output = [6.1734, 6.1710];
        %     testCase.verifyEqual(entropies{1}, expected_output(1), 'AbsTol', 0.05);
        %     testCase.verifyEqual(entropies{2}, expected_output(2), 'AbsTol', 0.05);         
        % end
        % 
        % function testInvalidX(testCase)
        %     rng('default');
        %     X = [1, 2, 3; 4, NaN, 6];
        %     Y = [1, 2, 3];
        %     opts.n_binsX = 10;
        %     opts.n_binsY = 5;
        %     opts.verbose = false;
        %     opts.bias = 'gsb';
        %     opts.method = 'gs';
        %     outputsList = {'HX', 'HXY'};  
        %     testCase.verifyError(@() entropy(X, Y, opts, outputsList), 'sanitycheck:NaNInput');
        % end
        % 
        % function testInvalidOutputList(testCase)
        %     nTrials = 100;
        %     nDimensionsX = 3;
        %     nDimensionsY = 1;
        %     X = randn(nDimensionsX, nTrials);
        %     Y = randi([1, 5], nDimensionsY, nTrials);
        %     opts.n_binsX = 10;
        %     opts.n_binsY = 5;
        %     opts.bias = 'gsb';
        %     opts.method = 'gs';
        %     outputsList = {'InvalidOutput'};
        %     testCase.verifyError(@() entropy(X, Y, opts, outputsList), 'entropy:InvalidOutputList');        
        % end  
        % 
        % function testInvalidBias(testCase)
        %     nTrials = 100;
        %     nDimensionsX = 3;
        %     nDimensionsY = 1;
        %     X = randn(nDimensionsX, nTrials);
        %     Y = randi([1, 5], nDimensionsY, nTrials);
        %     opts.n_binsX = 10;
        %     opts.n_binsY = 5;
        %     opts.verbose = false;
        %     opts.bias = 'Invalid';
        %     opts.method = 'gs';
        %     outputsList = {'HX', 'HXY'};
        %     testCase.verifyError(@() entropy(X, Y, opts, outputsList), 'entropy:UndefinedBiasCorrectionMethod');         
        % end 
        % 
        % function testInvalidBinning(testCase)
        %     nTrials = 100;
        %     nDimensionsX = 3;
        %     nDimensionsY = 1;
        %     X = randn(nDimensionsX, nTrials);
        %     Y = randi([1, 5], nDimensionsY, nTrials);
        %     opts.n_binsX = 10;
        %     opts.bin_methodY = 'eqpop';
        %     opts.bin_methodX = 'eqpop';
        %     opts.verbose = false;
        %     opts.bias = 'gsb';
        %     opts.method = 'gs';
        %     outputsList = {'HX', 'HXY'};
        %     testCase.verifyError(@() entropy(X, Y, opts, outputsList), 'entropy:UndefinedBinningY');      
        % end 
    end
end