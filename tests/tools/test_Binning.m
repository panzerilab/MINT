classdef test_Binning < matlab.unittest.TestCase


    methods (Test)
        function testEmptyInputGseqspace(testCase)
            functionCall = @() gseqspace();
            verifyError(testCase, functionCall, 'MATLAB:minrhs', 'Expected MATLAB:minrhs error.');
        end

        function testEmptyInputEqspace(testCase)
            functionCall = @() eqspace();
            verifyError(testCase, functionCall, 'MATLAB:minrhs', 'Expected MATLAB:minrhs error.');
        end

        function testEmptyInputEqpop(testCase)
            functionCall = @() eqpop();
            verifyError(testCase, functionCall, 'MATLAB:minrhs', 'Expected MATLAB:minrhs error.');
        end

        function testEmptyInputCeqspace(testCase)
            functionCall = @() ceqspace();
            verifyError(testCase, functionCall, 'MATLAB:minrhs', 'Expected MATLAB:minrhs error.');
        end

        function testGseqspace(testCase)
            rng('default');
            X = randn(100, 1);
            Nb = 10;
            edges = gseqspace(X, Nb);
            assert(length(edges) == Nb + 1, 'Number of edges is incorrect');
            assert(all(diff(edges) > 0), 'Edges are not in ascending order');

            X = randn(200, 1);
            Nb = 15;
            N = 4;
            opts.internalPar = N;
            edges = gseqspace(X, Nb, opts);
            assert(length(edges) == Nb + 1, 'Number of edges is incorrect');
            assert(all(diff(edges) > 0), 'Edges are not in ascending order');

            %Invalid Input
            X = randn(120, 1);
            Nb = 8;
            opts.internalPar = -2;
            try
                edges = gseqspace(X, Nb,opts); % Invalid N
                errorOccurred = false;
            catch ME
                errorOccurred = contains(ME.message, 'Parameter must be a positive scalar.');
            end
            testCase.verifyTrue(errorOccurred, 'Invalid N was not recognised');
        end

        function testEqspace(testCase)
            rng('default');
            X = randn(100, 1);
            nb = 10;
            edg = eqspace(X, nb);
            assert(length(edg) == nb + 1, 'Number of edges is incorrect');
            assert(all(diff(edg) > 0), 'Edges are not in ascending order');

            X = randn(200, 1);
            nb = 15;
            customRange = [-5, 5];
            opts.internalPar = customRange;
            edg = eqspace(X, nb, opts);
            assert(length(edg) == nb + 1, 'Number of edges is incorrect');
            assert(all(diff(edg) > 0), 'Edges are not in ascending order');
            assert(all(edg >= customRange(1)), 'Edges exceed the lower limit of the custom range');
            assert(all(edg <= customRange(2)+1), 'Edges exceed the upper limit of the custom range');

            %Invalid Input
            invalidRange = [2, -2];
            opts.internalPar = invalidRange;
            try
                edg = eqspace(X, nb,  opts);
                errorOccurred = false;
            catch ME
                errorOccurred = contains(ME.message, 'Invalid range.');
            end
            testCase.verifyTrue(errorOccurred, 'Invalid range was not recognised');
        end


        function testCeqspace(testCase)
            rng('default');
            X = randn(100, 1);
            nb = 10;
            edgs = ceqspace(X, nb);
            assert(length(edgs) == nb + 1, 'Number of edges is incorrect');
            assert(all(diff(edgs) > 0), 'Edges are not in ascending order');

            %Invalid Input
            X = randn(120, 1);
            nb = 8;
            invalidD = [2, -2]; % Invalid deviation
            opts.internalPar = invalidD;
            try
                edgs = ceqspace(X, nb, opts);
                errorOccurred = false;
            catch ME
                errorOccurred = contains(ME.message, 'Invalid parameter.');
            end
            testCase.verifyTrue(errorOccurred, 'Invalid deviation not recognised');
        end
       

        function testScottsRule(testCase)
            rng('default');
            data = randn(100, 1);
            num_bins = scottsRule(data);
            assert(num_bins == 8, 'Test failed: unexpected output');

            data = ones(50, 1);
            num_bins = scottsRule(data);
            assert(num_bins == 1, 'Test failed: unexpected output');
        end

        function testFreedmanDiaconisRule(testCase)
            rng('default');
            data = randn(100, 1);
            num_bins = freedmanDiaconisRule(data);
            assert(num_bins == 10, 'Test failed: unexpected output');

            data = ones(50, 1);
            num_bins = freedmanDiaconisRule(data);
            assert(num_bins == 1, 'Test failed: unexpected output');
        end
    end
end