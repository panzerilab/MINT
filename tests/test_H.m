classdef test_H < matlab.unittest.TestCase
 methods (Test)

        function testEmptyInput(testCase)
            testCase.verifyError(@()MI(), 'MI:notEnoughInput');
        end

        function testZerosInput(testCase)
            warning('off', 'all');
            S = zeros(1, 100);
            R = zeros(1, 100);
            H_result = H({R, S});
            assert(H_result{1} < 0.01, 'Value should be ~zero.');
            warning('on', 'all');
        end

        function testNaNInput(testCase)
            rng("default")
            S = zeros(1, 100);
            S(2, 5) = NaN;
            R = zeros(1, 100);           
            opts.NaN_handling = 'error';
            opts.suppressWarnings = true;
            testCase.verifyError(@()H({R, S}, opts), 'checkInputs:NaNDetected'); 
        end

        function testRandiValues(testCase)
            warning('off', 'all');
            rng("default")
            S = randi([1, 3], 1, 100);
            R = randi([1, 3], 1, 100);           
            H_result = H({R, S});
            assert((H_result{1} - log2(3) < 0.01) && (H_result{2} - log2(3) < 0.01), 'Value should be ~zero.');
            warning('on', 'all');
        end
       
    end

end