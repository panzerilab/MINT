classdef test_PID < matlab.unittest.TestCase

    methods (Test)
        function testEmptyInput(testCase)
            testCase.verifyError(@()PID(), 'PID:notEnoughInput');
        end

        function testZerosInput(testCase)
            warning('off', 'all');
            Y = zeros(3, 100);
            X1 = zeros(3, 100);
            X2 = zeros(3, 100);         
            PID_result = PID({X1; X2}, Y);
            assert(all(abs(PID_result) < 0.01), 'All values should be approximately zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            X1 = rand(3, 100);
            X2 = rand(3, 100);
            Y = rand(3, 100);
            Y(2, 5) = NaN;
            testCase.verifyError(@()PID({X1; X2}, Y), 'sanitycheck:NaNInput');
        end
      

        function testRandValues(testCase)
            rng('default');
            warning('off', 'all');
            Y = randi(2,1, 100);
            X1 = rand(1, 100);
            X2 = rand(1, 100);
            opts.bias = 'qe';
            opts.parallel = 1;
            opts.xtrp= 20;
            opts.bin_methodX = 'eqspace';
            opts.bin_methodY = 'eqpop';
            opts.n_binsX = 2;
            opts.n_binsY = 2;
            [PID_out, PID_uncorrected] = PID({X1; X2}, Y, opts);
            randidx = randi(100, 1, 100);
            [PID_outsh, PID_uncorrectedsh] = PID({X1; X2}, Y(:,randidx), opts);
            expected_output = [0. 0. 0. 0.]; %[1.6416, 1.2575, 1.2020, 0.5052];
            difference_within_tolerance = all(abs(PID_out - expected_output) < 0.1);
            assert(difference_within_tolerance, 'The output does not match the expected values within the specified tolerance.');
            warning('on', 'all');
        end


        function testValid_InvalidInputs(testCase)
            rng('default');
            warning('off', 'all');
            Y = randi(2, 1, 100) + .2;
            X1 = rand(1, 100);
            X2 = rand(1, 100);
            X = {X1; X2};
            opts.bias = 'qe';
            opts.parallel = 1;
            opts.xtrp= 20;
            opts.bin_methodX = 'eqpop';
            opts.bin_methodY = 'eqpop';
            opts.n_binsX = 2;
            opts.n_binsY = 2;
            opts.redundancy_measure = 'I_BROJA';
            opts.btsp=20;
            [PID_out, PID_uncorrected, ~, ~, ~, PID_btsp] = PID(X, Y, opts);
            expected_output = [0,0,0,0];

            difference_within_tolerance = all((PID_out - expected_output) < 0.1);
            assert(difference_within_tolerance, 'The output does not match the expected values within the specified tolerance.');


            opts.bias = 'invalid';
            testCase.verifyError(@()PID(X, Y, opts), 'PID:UndefinedBiasCorrectionMethod');

            opts.bias = 'naive';
            opts.bin_methodX  = 'invalid';
            testCase.verifyError(@()PID(X, Y, opts), 'PID:InvalidBinningMethod');
            warning('on', 'all');
        end
      
    end
end
