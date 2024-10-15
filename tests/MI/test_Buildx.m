classdef test_Buildx < matlab.unittest.TestCase

    methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()buildx(), 'buildx:notEnoughInput');
        end
        
        function testValidInput(testCase)
            rng('default');
            Y = randi([1, 5], 1, 100);
            X_in = randn(10, 100);
            [X, nt] = buildx(Y, X_in);
            testCase.verifyEqual(size(X), [10, max(nt), 5]);
            testCase.verifyEqual(size(nt), [5, 1]);
            for i = 1:5
                testCase.verifyEqual(sum(nt(i) == histc(Y, unique(Y))) > 0, true);
            end
        end

        function testYNotVectorError(testCase)
            Y = [1, 2; 3, 4]; % Not a vector, but a 2x2 matrix
            X_in = rand(5, 4); 
            functionCall = @() buildx(Y, X_in);
            verifyError(testCase, @() buildx(Y, X_in), 'buildx:StimNotAVector');
        end

        function testInvalidResponseArray(testCase)
            Y = [1, 2, 3];
            X_in = rand(4, 5); % 4 dimensions, but 5 trials
            verifyError(testCase, @() buildx(Y, X_in), 'buildx:differentTotNt');
        end
    end
end