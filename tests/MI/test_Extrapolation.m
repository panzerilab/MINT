classdef test_Extrapolation < matlab.unittest.TestCase

    methods (Test)

        % Helpers/xtrploop
        function testEmptyInputXtrploop(testCase)
            testCase.verifyError(@()xtrploop(), 'xtrploop:notEnoughInput');
        end
        
        function testRandInputXtrploop(testCase)
            rng('default');
            X = rand(3, 4, 5);          
            opts.nt = 4;
            opts.method = 'dr';
            opts.bias = 'naive';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);
            [HX, HXY, HlX, HlXY, HiX, HiXY, ChiX, HshX, HshXY] = xtrploop(X, pars);
            tolerance = 0.01;
            areAllZero = all([HX, HXY', HlX, HlXY, HiX, HiXY, ChiX, HshX, HshXY] <= tolerance);
            testCase.verifyTrue(areAllZero);
        end

        function testValidInputXtrploop(testCase)
            rng('default');
            numDimensions = 3;
            numTrials = 10;
            numStimuli = 3;
            X = zeros(numDimensions, numTrials, numStimuli);
            for i = 1:numStimuli
                meanVector = randn(numDimensions, 1) * i;
                covMatrix = eye(numDimensions);
                for j = 1:numTrials
                    X(:, j, i) = mvnrnd(meanVector, covMatrix)+ 5;
                end
            end
            opts.nt = 10;
            opts.method = 'dr';
            opts.bias = 'naive';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);
            [HX, HXY] = xtrploop(X, pars);
            notZero = all([HX, HXY'] > 0);
            testCase.verifyTrue(notZero);    
        end

        function testInvalidParsXtrploop(testCase)
            X = rand(3, 4, 5);
            pars = struct();
            errorOccurred = false;
            try
                [HX, HXY, HlX, HlXY, HiX, HiXY, ChiX, HshX, HshXY] = xtrploop(X, pars);  
            catch
                errorOccurred = true;
            end
            testCase.verifyTrue(errorOccurred);
        end

        % Quadratic/qadratic_extrapolation
        function testEmptyInputQE(testCase)
            testCase.verifyError(@()quadratic_extrapolation(), 'quadraticExtrpolation:notEnoughInput');
        end

        function testRandInputQE(testCase)
        rng('default');
            X = rand(3, 4, 5);          
            opts.nt = 4;
            opts.method = 'dr';
            opts.bias = 'naive';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);
            [HX, HXY, HlX, HlXY, HiX, HiXY, ChiX, HshX, HshXY] = quadratic_extrapolation(X, pars);
            tolerance = 0.01;
            areAllZero = all([HX, HXY', HlX, HlXY, HiX, HiXY, ChiX, HshX, HshXY] <= tolerance);
            testCase.verifyTrue(areAllZero);
        end

        function testValidInputQE(testCase)
            rng('default');
            numDimensions = 3;
            numTrials = 10;
            numStimuli = 3;
            X = zeros(numDimensions, numTrials, numStimuli);
            for i = 1:numStimuli
                meanVector = randn(numDimensions, 1) * i;
                covMatrix = eye(numDimensions);
                for j = 1:numTrials
                    X(:, j, i) = mvnrnd(meanVector, covMatrix)+ 5;
                end
            end
            opts.nt = 10;
            opts.method = 'dr';
            opts.bias = 'naive';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);
            [HX, HXY] = quadratic_extrapolation(X, pars);
            notZero = all([HX, HXY'] > 0);
            testCase.verifyTrue(notZero);    
        end
    end
end