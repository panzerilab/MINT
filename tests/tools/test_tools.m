classdef test_tools < matlab.unittest.TestCase

    methods (Test)     
        function reduceDim(testCase)
            testCase.verifyError(@()reduce_dim(), 'reduce_dim:notEnoughInput');

            %Valid Input
            Nd_1d = [1, 2, 3, 4];
            expected_output_1d = Nd_1d;
            assert(isequal(reduce_dim(Nd_1d, 1), expected_output_1d),'Test failed: unexpected output');

            Nd_1d = [1, 2, 1; 2, 3, 1];
            expected_output_1d = [3, 6, 1];
            assert(isequal(reduce_dim(Nd_1d, 1), expected_output_1d),'Test failed: unexpected output');
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
            testCase.verifyError(@()shuffle_core(), 'shuffle_core:notEnoughInput');
            
            % Valid Input
            rng('default');
            behav_data = randi([1, 10], 100, 1);
            neural_data = rand(100, 10); 
            consistency = 0;
            index = [1, 0];
            shuffled_data = shuffle_core(behav_data, neural_data, consistency, index);
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
            shuffled_data = shuffle_core(behav_data, neural_data, consistency, index);
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
            shuffle_core(behav_data, neural_data, consistency, index);
            [warnMsg, warnId] = lastwarn;
            expectedWarningMsg = 'Shuffling the second dimension (neuralObjects) may alter probability distributions of individual objects.';
            assert(contains(warnMsg, expectedWarningMsg), 'Expected warning message was not issued.');
        end
    end
end