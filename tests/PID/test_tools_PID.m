classdef test_tools_PID < matlab.unittest.TestCase

    methods (Test)

        function testBuildP(testCase)
            rng('default');
            warning('off', 'all');
            % empty input
            testCase.verifyError(@()build_p(), 'buildp:notEnoughInput');

            % zeros input
            S = zeros(1, 100);
            R = zeros(1, 100);
            C = zeros(1, 100);
            [p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] = build_p(S, R, C);
            expected = [1, 1, 1, 1, 1, 3];
            difference = all(abs([p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] - expected) == 0);
            assert(difference);

            % NaN Input
            S(1,1) = NaN;
            testCase.verifyError(@()build_p(S, R, C), 'buildp:NaNInput');

            % Valid Input
            S = [1; 2; 1; 2];
            R = [3; 3; 4; 4];
            C = [5; 6; 5; 6];
            [p_src, p_crs, n_S, n_R, n_C, n_singleton_dims] = build_p(S, R, C);
            assert(isequal(size(p_src), [2, 2, 2]));
            assert(isequal(n_S, 2) && isequal(n_R, 2) && isequal(n_C, 2));
            warning('on', 'all');
        end

        function testCalculatePid(testCase)
            warning('off', 'all');
            rng('default');
            % empty input
            testCase.verifyError(@()calculate_pid(), 'calculatepid:notEnoughInput');

            % zeros input
            Y = zeros(3, 100);
            X1 = zeros(3, 100);
            X2 = zeros(3, 100);
            n_trials = 100;
            n_split_trials = 10;
            PID = calculate_pid(Y,X1,X2,n_trials,n_split_trials);
            assert(all(abs(PID) < 0.01), 'All values should be approximately zero.');

            % NaN Input
            Y(1,1) = NaN;
            testCase.verifyError(@()calculate_pid(Y,X1,X2,n_trials,n_split_trials), 'calculatepid:NaNInput');

            % Valid Input
            X1 = rand(3, 100);
            X2 = rand(3, 100);
            Y = rand(3, 100);
            n_trials = 100;
            n_split_trials = 10;
            PID = calculate_pid(Y,X1,X2,n_trials,n_split_trials);
            expected_output = [3.3219, 0, 0, 0];
            difference_within_tolerance = all(abs(PID - expected_output) < 0.01);
            assert(difference_within_tolerance, 'The output does not match the expected values within the specified tolerance.');
            warning('on', 'all');
        end

        function testCreateProbTs(testCase)
            warning('off', 'all');
            rng('default');
            % empty input
            testCase.verifyError(@()create_prob_ts(), 'createprobts:notEnoughInput');

            % zeros input
            p_distr = 0;
            sources = 1;
            p_ts = create_prob_ts(p_distr, sources);
            assert(p_ts == 0, 'p_ts should be zero.');

            % NaN Input
            p_distr = NaN;
            testCase.verifyError(@()create_prob_ts(p_distr, sources), 'createprobts:NaNInput');

            % Valid Input
            p_distr = rand(2,3,4);
            p_distr = p_distr/sum(p_distr(:));
            sources = 2;
            p_ts = create_prob_ts(p_distr, sources);
            warning('on', 'all');
        end

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

        function testPidExtrapolation(testCase)
            warning('off', 'all');
            rng('default');
            %empty input
            testCase.verifyError(@()extrapolation(), 'extrapolation:notEnoughInput');

            % zeros input
            Y_b = ones(3, 100);
            X1_b = ones(3, 100);
            X2_b = ones(3, 100);
            opts.xtrp = 4;
            opts.function = @pidimin;
            opts.parallel = 0;
            opts.bias = 'qe';
            btspflag = false;
            opts.n_trials = 100;
            [PID_corrected, PID_naive] = extrapolation(@pid_core, [Y_b; X1_b; X2_b], opts, btspflag);
            %[PID_corrected, PID_naive] = extrapolation(Y_b, X1_b, X2_b, opts, btspflag);
            equals_zero = all(PID_corrected == 0) & all(PID_naive == 0);
            assert(equals_zero);

            % NaN Input
            Y_b(1,1) = NaN;
            testCase.verifyError(@()extrapolation(@pid_core, [Y_b; X1_b; X2_b], opts, btspflag), 'sanitycheck:NaNInput');

            % Valid Input
            Y_b = randi([1,2], 1, 100);
            X1_b = randi([1,2], 1, 100);
            X2_b =randi([1,2], 1, 100);
            [PID_corrected, PID_naive] = extrapolation(@pid_core, [Y_b; X1_b; X2_b], opts, btspflag);
            equals_zero = all(abs(PID_corrected) < 0.1) & all(PID_naive < 0.1);
            assert(equals_zero);
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


