classdef test_tools < matlab.unittest.TestCase

    methods (Test)

      function testCheckOpts(testCase)
            testCase.verifyError(@()check_opts(), 'checkopts:notEnoughInput');

            %Valid Input
            opts = struct('bin_methodX1', 'eqspace', 'n_binsX1', 10, ...
                'bias', 'naive', 'max_draws_per_split_number', 5, ...
                'btsp', 0, 'nsources', 2);
            measure_name = 'PID';
            testCase.verifyWarningFree(@()check_opts(opts, measure_name));

            %Invalid Input
            opts = struct('bin_methodX1', 'invalid', 'n_binsX1', 10, ...
                'bias', 'naive', 'max_draws_per_split_number', 5, ...
                'btsp', 0, 'nsources', 2);
            measure_name = 'PID';
            testCase.verifyError(@()check_opts(opts, measure_name), 'checkOpts:invalidBinMethod')

            opts = struct('bin_methodX1', 'eqspace', 'n_binsX1', 10, ...
                'bias', 'naive', 'max_draws_per_split_number', 5, ...
                'btsp', 0, 'nsources', 3, 'redundancy_measure', 'I_BROJA');
            measure_name = 'PID';
            testCase.verifyError(@()check_opts(opts, measure_name), 'checkOpts:InvalidRedundancyMeasure')
        end

        function testConsitencyCheck(testCase)
            testCase.verifyError(@()consistency_check(), 'consistencycheck:notEnoughInput');

            %Valid Input
            X1 = rand(5, 100, 100);
            X2 = rand(5, 100, 100);
            ntrials = consistency_check({X1, X2});
            assert(ntrials == 100);

            %Invalid Input
            X1 = rand(5, 100, 100);
            X2 = rand(5, 90, 10);
            testCase.verifyError(@()consistency_check({X1, X2}), 'consistencycheck:InvalidInput');
        end

        function testMapNdArrayTo1D(testCase)
            testCase.verifyError(@()map_Nd_array_to_1d(), 'mapNdarrayto1d:notEnoughInput');

            %Valid Input
            Nd_1d = [1, 2, 3, 4];
            expected_output_1d = Nd_1d;
            assert(isequal(map_Nd_array_to_1d(Nd_1d), expected_output_1d),'Test failed: unexpected output');

            Nd_1d = [1, 2, 1; 2, 3, 1];
            expected_output_1d = [3, 6, 1];
            assert(isequal(map_Nd_array_to_1d(Nd_1d), expected_output_1d),'Test failed: unexpected output');

            %Invalid Input
            Nd_1d = [1, 2, 1; 1, 2, NaN];
            testCase.verifyError(@()map_Nd_array_to_1d(Nd_1d), 'mapNdarrayto1d:NaNInput');
        end

        function testMarginalizePdf(testCase)
            rng('default');
            testCase.verifyError(@()marginalize_pdf(), 'marginalizedpdf:notEnoughInput');

            %Valid Input
            pdf_size = [10, 20, 15];
            pdf = rand(pdf_size);
            survivingdims_single = 2;
            survivingdims_vector = [1, 3];
            iscont_vector = logical([1, 0, 1]);
            xi_1 = linspace(-2, 2, pdf_size(1));
            xi_2 = linspace(-3, 3, pdf_size(2));
            xi_3 = linspace(-1, 1, pdf_size(3));
            marginalized_pdf_single = marginalize_pdf(pdf, survivingdims_single, iscont_vector, xi_1, xi_2, xi_3);
            marginalized_pdf_vector = marginalize_pdf(pdf, survivingdims_vector, iscont_vector, xi_1, xi_2, xi_3);
            assert(all(size(marginalized_pdf_single) == [pdf_size(survivingdims_single), length(survivingdims_single)]), 'Test failed: Invalid size for marginalized_pdf_single');
            assert(all(size(marginalized_pdf_vector) == [pdf_size(survivingdims_vector)]), 'Test failed: Invalid size for marginalized_pdf_single');

            %Invalid Input
            pdf_size = [10, 20, 15];
            pdf = rand(pdf_size);
            survivingdims = [1, 1];
            iscont_vector = logical([1, 0, 1]);
            xi_1 = linspace(-2, 2, pdf_size(1));
            xi_2 = linspace(-3, 3, pdf_size(2));
            xi_3 = linspace(-1, 1, pdf_size(3));
            testCase.verifyError(@()marginalize_pdf(pdf, survivingdims, iscont_vector, xi_1, xi_2, xi_3), 'marginalizedpdf:RepeatedIndices');

            survivingdims = [1, 3];
            testCase.verifyError(@()marginalize_pdf(pdf, survivingdims, iscont_vector, xi_1, xi_2), 'marginalizedpdf:MismatchArguments');

            pdf(1,3) = NaN;
            testCase.verifyError(@()marginalize_pdf(pdf, survivingdims, iscont_vector, xi_1, xi_2, xi_3), 'marginalizedpdf:NaNInput');
        end


        function testNonPoissonSpikeGen(testCase)
            testCase.verifyError(@()non_poisson_spike_gen(), 'nonpoissonspikegen:notEnoughInput');

            %Valid Input
            time = 0:0.1:12;
            rate = 6;
            noise_prob = 0;
            spikes = non_poisson_spike_gen(time, rate, noise_prob);
            assert(length(spikes) == length(time), 'Test failed: Incorrect length of spikes');
            assert(all(ismember(spikes, [0, 1])), 'Test failed: spikes should contain only 0s and 1s');

            time = 0:0.1:10;
            rate = 8;
            noise_prob = 0.2;
            spikes = non_poisson_spike_gen(time, rate, noise_prob);
            assert(length(spikes) == length(time), 'Test failed: Incorrect length of spikes');
            assert(all(ismember(spikes, [0, 1])), 'Test failed: spikes should contain only 0s and 1s');

            %Invalid Input
            time = 0;
            rate = 6;
            noise_prob = 0;
            testCase.verifyError(@()non_poisson_spike_gen(time, rate, noise_prob), 'MATLAB:badsubscript');

            time = NaN;
            testCase.verifyError(@()non_poisson_spike_gen(time, rate, noise_prob), 'nonpoissonspikegen:NaNInput');
        end

        % function testPdf(testCase)
        %     rng('default');
        %     testCase.verifyError(@()pdf(), 'pdf:notEnoughInput');
        %
        %     data = rand(3, 50, 2);
        %     opts.method = 'none';
        %     probdist = pdf(data, opts);
        %
        %
        %     data = rand(3, 50, 2);
        %     opts.nb = 5;
        %     opts.method = 'eqspace';
        %     probdist = pdf(data, opts);
        %     assert(isequal(size(probdist), [5, 5, 5]), 'Test failed: Incorrect Output Dimension');
        %     assert(isequal(sum(sum(sum(probdist))), 1),'Test failed: Sum of all probabilities must be 1');
        %
        %     %Invalid Input
        %     data = rand(4, 30, 3);
        %     opts_invalid.nb = 'invalid';
        %     try
        %         probdist = pdf(data, opts_invalid);
        %         error('Test Case failed: Invalid opts.nb was not detected.');
        %     catch
        %     end
        % end

        function testpoissonSpikeGen(testCase)
            testCase.verifyError(@()poisson_spike_gen(), 'poissonspikegen:notEnoughInput');

            %Valid Input
            time = 0:0.1:12;
            rate = 6;
            noise_prob = 0;
            spikes = poisson_spike_gen(time, rate, noise_prob);
            assert(length(spikes) == length(time), 'Test failed: Incorrect length of spikes');
            assert(all(ismember(spikes, [0, 1])), 'Test failed: spikes should contain only 0s and 1s');

            time = 0:0.1:10;
            rate = sin(time);
            noise_prob = 0;
            spikes = poisson_spike_gen(time, rate, noise_prob);
            assert(length(spikes) == length(time), 'Test failed: Incorrect length of spikes');
            assert(all(ismember(spikes, [0, 1])), 'Test failed: spikes should contain only 0s and 1s');

            time = 0:0.1:10;
            rate = 8;
            noise_prob = 0.2;
            spikes = poisson_spike_gen(time, rate, noise_prob);
            assert(length(spikes) == length(time), 'Test failed: Incorrect length of spikes');
            assert(all(ismember(spikes, [0, 1])), 'Test failed: spikes should contain only 0s and 1s');

            %Invalid input
            time = 0;
            rate = 6;
            noise_prob = 0;
            testCase.verifyError(@()poisson_spike_gen(time, rate, noise_prob), 'MATLAB:badsubscript');

            time = NaN;
            testCase.verifyError(@()poisson_spike_gen(time, rate, noise_prob), 'poissonspikegen:NaNInput');
        end


        function testShuffle(testCase)
            testCase.verifyError(@()shuffle(), 'shuffle:notEnoughInput');
            
            % Valid Input
            rng('default');
            behav_data = randi([1, 10], 100, 1);
            neural_data = rand(100, 10); 
            consistency = 0;
            index = [1, 0];
            shuffled_data = shuffle(behav_data, neural_data, consistency, index);
            assert(~isequal(neural_data, shuffled_data), 'Data should be shuffled');
            nbins = 10;
            [counts_neural, ~] = histcounts(neural_data(:,1), nbins);
            [counts_shuffled, ~] = histcounts(shuffled_data(:,1), nbins);
            assert(isequal(counts_neural, counts_shuffled));
            uniqueBehav = unique(behav_data);
            for i = 1:length(uniqueBehav)
                currentBehav = uniqueBehav(i); 
                idx = behav_data == currentBehav; 
                currentNeuralData = neural_data(idx, :);
                currentShuffledData = shuffled_data(idx, :);
                for object_idx = 1:size(neural_data, 2)
                    [counts_neural, ~] = histcounts(currentNeuralData(:,object_idx), nbins);
                    [counts_shuffled, ~] = histcounts(currentShuffledData(:,object_idx), nbins);
                    assert(~isequal(counts_neural, counts_shuffled), 'Histograms of the first dimension should differ.');
                end
            end

            consistency = 1;
            shuffled_data = shuffle(behav_data, neural_data, consistency, index);
            assert(~isequal(neural_data, shuffled_data), 'Data should be shuffled');
            nbins = 10;
            [counts_neural, ~] = histcounts(neural_data(:,1), nbins);
            [counts_shuffled, ~] = histcounts(shuffled_data(:,1), nbins);
            assert(isequal(counts_neural, counts_shuffled));
            uniqueBehav = unique(behav_data);
            for i = 1:length(uniqueBehav)
                currentBehav = uniqueBehav(i); 
                idx = behav_data == currentBehav; 
                currentNeuralData = neural_data(idx, :);
                currentShuffledData = shuffled_data(idx, :);
                for object_idx = 1:size(neural_data, 2)
                    [counts_neural, ~] = histcounts(currentNeuralData(:,object_idx), nbins);
                    [counts_shuffled, ~] = histcounts(currentShuffledData(:,object_idx), nbins);
                    assert(isequal(counts_neural, counts_shuffled), 'Shuffling should be consistent across behavData.');
                end
            end

            % Test Warning
            lastwarn('');
            disp('Dont worry, this warning is part of the test ;)');
            behav_data = randi([1, 10], 100);
            neural_data = rand(100, 10); 
            consistency = 0;
            index = [1, 1];
            shuffle(behav_data, neural_data, consistency, index);
            [warnMsg, warnId] = lastwarn;
            expectedWarningMsg = 'Shuffling the second dimension (neuralObjects) may alter probability distributions of individual objects.';
            assert(contains(warnMsg, expectedWarningMsg), 'Expected warning message was not issued.');
        end



        function testWindowing(testCase)
            testCase.verifyError(@()windowing(), 'windowing:notEnoughInput');

            %Valid Inputs;
            time = 1:10;
            nWindows = 2;
            window = windowing(time, nWindows);
            assert(isequal(window, [5, 5]), 'Test failed');

            time = [1, 2, 4, 6, 9, 10];
            nWindows = 3;
            window = windowing(time, nWindows);
            assert(isequal(window, [2, 2, 2]), 'Test failed');

            time = 1:10;
            window = windowing(time, nWindows);
            assert(isequal(window, [3, 3, 4]), 'Test failed');

            time = 0.1:0.1:1;
            window = windowing(time, nWindows);
            assert(isequal(window, [3, 3, 4]), 'Test failed');

            %Invalid Input
            time = 1:5;
            nWindows = 10;            
            testCase.verifyError(@()windowing(time, nWindows), 'windowing:InvalidNumberOfWindows');

            time = 1:10;
            nWindows = 0;
            testCase.verifyError(@()windowing(time, nWindows), 'windowing:InvalidNumberOfWindows');

            nWindows = NaN;
            testCase.verifyError(@()windowing(time, nWindows), 'windowing:InvalidNumberOfWindows');
        end
    end
end