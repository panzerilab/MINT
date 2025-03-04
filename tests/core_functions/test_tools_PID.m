classdef test_tools_PID < matlab.unittest.TestCase

    methods (Test)

        function testIIBroja(testCase)
            warning('off', 'all');
            rng('default');
            % empty input
            testCase.verifyError(@()iiBROJA(), 'iibroja:notEnoughInput');

            % zero input
            pdf = 0;
            testCase.verifyError(@()iiBROJA(pdf), 'iibroja:InvalidInput');

            % negative input
            pdf = -1;
            testCase.verifyError(@()iiBROJA(pdf), 'iibroja:InvalidInput');

            % NaN Input
            pdf = NaN;
            testCase.verifyError(@()iiBROJA(pdf), 'iibroja:NaNInput');

            % Valid Input
            pdf = rand(3, 3, 3);
            ii_result = iiBROJA(pdf);
            expected_ii_result = [0.0082, 0.0082, 0.0459, 0, 0.1178, 0.0158, 0.0161,0.0382, 0.1254];
            assert(length(ii_result) == 9);
            assert(abs(sum(expected_ii_result)-sum(ii_result)) < 0.01);
            warning('on', 'all');
        end

        function testPIDBroja(testCase)
            warning('off', 'all');
            % empty input
            testCase.verifyError(@()pidBROJA(), 'pidbroja:notEnoughInput');

            % zero input
            pdf = 0;
            testCase.verifyError(@()pidBROJA(pdf), 'pidbroja:InvalidInput');

            % negative input
            pdf = -1;
            testCase.verifyError(@()pidBROJA(pdf), 'pidbroja:InvalidInput');

            % NaN Input
            pdf = NaN;
            testCase.verifyError(@()pidBROJA(pdf), 'pidbroja:NaNInput');

            % Valid Input
            pdf = zeros(2, 2, 2);
            pdf(1, 1, 1) = 0.5;
            pdf(2, 2, 2) = 0.5;
            pid_v = pidBROJA(pdf);
            expected_pid_v = [1, 0, 0, 0];
            res = sum(round(pid_v - expected_pid_v));           
            verifyEqual(testCase, res, 0, 'AbsTol', 1e-4, 'RelTol', 1e-4);
            warning('on', 'all');
        end

        function testPidImin(testCase)
            warning('off', 'all');
            rng('default');
            % empty input
            testCase.verifyError(@()pidimin(), 'pidimin:notEnoughInput');

            % zero input
            pdf = 0;
            testCase.verifyError(@()pidimin(pdf), 'pidimin:InvalidInput');

            % negative input
            pdf = -1;
            testCase.verifyError(@()pidimin(pdf), 'pidimin:InvalidInput');

            % NaN Input
            pdf = NaN;
            testCase.verifyError(@()pidimin(pdf), 'pidimin:NaNInput');

            % Valid Input
            pdf = rand(3, 3, 3);
            pid_v = pidimin(pdf);
            expected_result = [0.0248, 0.0071, 0.0292, 0.1345];
            assert(length(pid_v) == 4);
            assert(abs(sum(expected_result)-sum(pid_v)) < 0.01);
            warning('on', 'all');
        end
    end
end


