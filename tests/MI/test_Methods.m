classdef test_Methods < matlab.unittest.TestCase

    methods (Test)
        
        function testEmptyInputXtrploop(testCase)
            testCase.verifyError(@()gaussian_method(), 'gaussianmethod:notEnoughInput');
        end
        
        function testGaussianMethodValidInput(testCase)
            rng('default');
            X = rand(3, 4, 5);          
            opts.nt = 4;
            opts.method = 'gs';
            opts.bias = 'gsb';
            responseMatrixName = 'responseMatrix';
            outputsList = {'hx', 'hxy'};
            opts.verbose = false;
            pars = build_parameters_structure(X, opts, responseMatrixName, outputsList);       
            [HX, HXY, HlX, HlXY, HiX, HiXY, ChiX, HshX, HshXY] = gaussian_method(X, pars);
            assert(all(isnan([ChiX, HiX, HiXY])));
            assert(all([HlX, HlXY, HshX, HshXY]) == 0);
            assert(all([HX, HXY']) > 0);
        end
    end
end