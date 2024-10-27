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
            PID_result = PID({X1 X2, Y});
            assert(all(cell2mat(PID_result) < 0.01), 'All values should be approximately zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            X1 = rand(3, 100);
            X2 = rand(3, 100);
            Y = rand(3, 100);
            Y(2, 5) = NaN;
            opts.NaN_handling = 'error';
            opts.supressWarnings = true;
            testCase.verifyError(@()PID({X1 X2, Y}, opts), 'checkInputs:NaNDetected'); 
        end
      
        function testRandValues(testCase)
            rng('default');            
            warning('off', 'all');
            Y = randi(2,1, 100);
            X1 = rand(1, 100);
            X2 = rand(1, 100);
            opts.bias = 'qe';
            opts.xtrp= 20;
            opts.supressWarnings = true;
            opts.bin_method = {'eqspace', 'eqpop', 'none'};
            opts.n_bins = {2, 2};
            PID_out = PID({X1, X2, Y}, opts); 
            assert(all(cell2mat(PID_out) < 0.1), 'All values should be approximately zero.');
            warning('on', 'all');
        end

        function testValid_InvalidInputs(testCase)
            rng('default');
            warning('off', 'all');
            Y = randi(2,1, 100);
            X1 = rand(1, 100);
            X2 = rand(1, 100);
            opts.bias = 'invalid';
            opts.xtrp= 20;
            opts.supressWarnings = true;
            opts.bin_method = {'eqspace', 'eqpop', 'none'};
            opts.n_bins = {2, 2};
            testCase.verifyError(@()PID({X1, X2, Y}, opts), 'Correction:UndefinedFunction');

            opts.bias = 'naive';
            opts.bin_method  = {'invalid'};
            testCase.verifyError(@()PID({X1, X2, Y}, opts), 'Binning:UndefinedFunction');
            warning('on', 'all');
        end
        
          function test_qdistrs(testCase)
            rng('default');
            warning('off', 'all');
            Y = randi(2, 1, 100);
            X1 = rand(1, 100);
            X2 = rand(1, 100);
            outputs = {'q_dist'};
            opts.redundancy_measure = 'I_BROJA';
            opts.bias = 'naive';            
            opts.supressWarnings = true;
            opts.bin_method = {'eqspace', 'eqpop', 'none'};
            opts.n_bins = {2, 2};
            PID_out = PID({X1, X2, Y}, outputs, opts); 
            warning('on', 'all');
        end

    end
end
