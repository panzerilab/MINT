classdef test_InputParsing < matlab.unittest.TestCase

    methods (Test)
        % BuildParameterStructure
        function testEmptyInputBuildParameterStructure(testCase)
            testCase.verifyError(@()build_parameters_structure(), 'buildParametersStructure:notEnoughInput');
        end

        function testValidInputBuildParameterStructure(testCase)
            rng('default');          
            X = rand(3, 4, 5);          
            opts.nt = 4;
            opts.method = 'dr';
            opts.bias = 'naive';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);
       
            verifyEqual(testCase, size(pars.Nt), [5, 1]);
            verifyEqual(testCase, pars.methodNum, 1);
            verifyEqual(testCase, pars.biasCorrNum, 0);
            verifyEqual(testCase, pars.numberOfSpecifiedOptions, 0);
            verifyEqual(testCase, pars.doHX, true);
            verifyEqual(testCase, pars.doHXY, true);
            verifyEqual(testCase, pars.addChecks, false);
        end

        function testInvalidNtBuildParameterStructure(testCase)
            rng('default');
            X = rand(2, 3, 4);
            opts.method = 'dr';
            opts.bias = 'naive';
            opts.nt = 'invalid';
            opts.verbose = false;
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            testCase.verifyError(@()build_parameters_structure(X, opts, responseMatrixName, outputsList), 'buildParametersStructure:InvalidInput');

            opts.nt = 5;
            testCase.verifyError(@()build_parameters_structure(X, opts, responseMatrixName, outputsList), 'buildParametersStructure:InvalidInput');
        end

        function testInvalidOutputBuildParameterStructure(testCase)
            rng('default');
            X = rand(2, 5, 4);
            opts.method = 'dr';
            opts.bias = 'naive';
            opts.nt = 5;
            opts.verbose = false;
            responseMatrixName = 'responseMatrix';
            outputsList = {};
            testCase.verifyError(@()build_parameters_structure(X, opts, responseMatrixName, outputsList), 'buildParametersStructure:InvalidOutputList');
        end

        function testInvalidMethodBuildParameterStructure(testCase)
            rng('default');
            X = rand(2, 3, 4);
            opts.method = 'invalid';
            opts.bias = 'naive';
            opts.nt = 3;
            opts.verbose = false;
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            testCase.verifyError(@()build_parameters_structure(X, opts, responseMatrixName, outputsList), 'buildParametersStructure:InvalidInput');
        end

        % AdditionalChecks
        function testEmptyInputAdditionalChecks(testCase)
            testCase.verifyError(@()additional_checks(), 'additionalchecks:notEnoughInput');
        end

        function testValidInputAdditionalChecks(testCase)
            warning('off', 'all');
            X = randi([0, 9], 3, 4, 5);  
            opts.nt = 4;
            opts.method = 'dr';
            opts.bias = 'naive';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);
            additional_checks(X, pars, opts);
            warning('on', 'all');
        end

        function testInvalidInputAdditionalChecks(testCase)
            X = rand(3, 4, 5);
            opts.nt = 4;
            opts.method = 'dr';
            opts.bias = 'naive';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);
            testCase.verifyError(@()additional_checks(X, pars, opts), 'additionalchecks:InvalidInput');            
        end
    end
end