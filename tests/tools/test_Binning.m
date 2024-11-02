classdef test_Binning < matlab.unittest.TestCase


    methods (Test)
       
          function testEqpop(testCase)
            rng('default');
            X = randi([1, 10],1, 100);
            Y = rand(2, 20, 100);
            opts.n_bins = {5, [6, 2]};
            opts.bin_method = {'eqpop'};
            binned_data = binning({X, Y}, opts);
            uniqueX = unique(binned_data{1});
            uniqueY_1 = unique(binned_data{1, 2}(1,:,:));
            uniqueY_2 = unique(binned_data{1, 2}(2,:,:));
            assert(length(unique(uniqueX)) == 5, 'Number of bins is incorrect');
            assert(length(unique(uniqueY_1)) == 6, 'Number of bins is incorrect');
            assert(length(unique(uniqueY_2)) == 2, 'Number of bins is incorrect'); 

            X = rand(1, 200);
            opts.bin_method = {'none'};
            binned_data = binning({X}, opts);
        end


        function testEqspace(testCase)
           rng('default');
            X = randi([1, 10],1, 100);
            Y = rand(2, 20, 100);
            opts.n_bins = {5, [6, 2]};
            opts.bin_method = {'eqspace'};
            binned_data = binning({X, Y}, opts);
            uniqueX = unique(binned_data{1});
            uniqueY_1 = unique(binned_data{1, 2}(1,:,:));
            uniqueY_2 = unique(binned_data{1, 2}(2,:,:));
            assert(length(unique(uniqueX)) == 5, 'Number of bins is incorrect');
            assert(length(unique(uniqueY_1)) == 6, 'Number of bins is incorrect');
            assert(length(unique(uniqueY_2)) == 2, 'Number of bins is incorrect'); 
        end

         function testUserEdges(testCase)
           rng('default');
            X = randi([1, 10],1, 100);
            Y = rand(2, 20, 100);
            opts.edges = {[1,3,5,7], [0.2, 0.6]};
            opts.bin_method = {'userEdges'};
            binned_data = binning({X, Y}, opts);
            uniqueX = unique(binned_data{1});
            uniqueY_1 = unique(binned_data{1, 2}(1,:,:));
            uniqueY_2 = unique(binned_data{1, 2}(2,:,:));
            assert(length(unique(uniqueX)) == 4, 'Number of bins is incorrect');
            assert(length(unique(uniqueY_1)) == 3, 'Number of bins is incorrect');
            assert(length(unique(uniqueY_2)) == 3, 'Number of bins is incorrect'); 
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