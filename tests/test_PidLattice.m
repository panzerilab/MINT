classdef test_PidLattice < matlab.unittest.TestCase

    methods (Test)

        function testConstructor(testCase)
            load('lat3sources.mat');
            nsources = 3;
            obj = pid_lattice(nsources);
            testCase.verifyEqual(obj.lat, lat3sources);
        end

        function testPowerSet(testCase)
            obj = pid_lattice(3); 
            var_list = [1, 2, 3];
            p_set = obj.power_set(var_list);
            expectedNumSubsets = 2^numel(var_list) - 1;
            testCase.verifyEqual(numel(p_set), expectedNumSubsets);
        end

        function testIsNodeRed(testCase)
            obj = pid_lattice(3);
            x = {[1,3],3};
            testCase.verifyFalse(is_node_red(obj, x));

            x = {[1,2],3};
            testCase.verifyTrue(is_node_red(obj, x));

            obj = pid_lattice(2);
            x = {[1,3],3};
            testCase.verifyFalse(is_node_red(obj, x));

            x = {[1,2]};
            testCase.verifyTrue(is_node_red(obj, x));

            x = {[1,3]};
            testCase.verifyFalse(is_node_red(obj, x));
        end

        function testNodeIsSubsetAny(testCase)
            obj = pid_lattice(3);
            x = {[1,3],3};
            i = 2;
            j = 1;
            testCase.verifyFalse(node_issubsetany(obj, x, i, j));
            
            x = {[1,2],3};
            testCase.verifyTrue(node_issubsetany(obj, x, i, j));
        end 

         function testNodeIsSubset(testCase)
            obj = pid_lattice(3);
            x = {1, 2};
            y = {1, 2};
            testCase.verifyTrue(node_issubset(obj, x, y));
            z = {4};
            testCase.verifyFalse(node_issubset(obj, x, z));
         end

         function testGetDown(testCase)
             obj = pid_lattice(3);
             node = {1,2,3};
             down_list = get_down(obj, node);
             testCase.verifyEqual(down_list{1, 1}{1, 1}, 1);
             testCase.verifyEqual(down_list{1, 1}{1, 2}, 2);
             testCase.verifyEqual(down_list{1, 1}{1, 3}, 3);

             obj = pid_lattice(2);
             node = {1,2,3};
             down_list = get_down(obj, node);
             testCase.verifyTrue(isempty(down_list));
         end

         function testGetStrictDown(testCase)
             obj = pid_lattice(3);
             node = {[1,2],3};
             down_list = get_strict_down(obj, node);
             testCase.verifyEqual(down_list{1, 1}{1, 1}, 1);
             testCase.verifyEqual(down_list{1, 1}{1, 2}, 3);
             testCase.verifyEqual(down_list{1, 2}{1, 1}, 2);

             node = {1,2,3};
             down_list = get_strict_down(obj, node);
             testCase.verifyTrue(isempty(down_list));

             obj = pid_lattice(2);
             node = {1,2,3};
             down_list = get_strict_down(obj, node);
             testCase.verifyTrue(isempty(down_list));
         end

         function testNodeIsSame(testCase)
             obj = pid_lattice(3);
             x = {1, 2};
             y = {1, 2};
             testCase.verifyTrue(node_issame(obj, x, y));
             
             y = {1,2,3};
             testCase.verifyFalse(node_issame(obj, x, y));             
         end

         % function testGetRed(testCase)
         %     obj = pid_lattice(3);
         %     get_red(obj);  
         % end
         % 
         % function testGetPi(testCase)
         %     obj = pid_lattice(3);
         %     pdf = rand(4,100,7);
         %     get_pi(obj, pdf)
         % end

         % function testImin(testCase)
         %     
         % end

        function testSpecificInformation(testCase)
            rng("default")
            obj = pid_lattice(3); 
            p_distr = [0.25, 0.25; 0.25, 0.25];
            spec_val_dim = 1;
            spec_val_index = 1;
            spec_info = obj.specific_information(p_distr, spec_val_dim, spec_val_index);
            expectedSpecInfo = 0; 
            testCase.verifyEqual(spec_info, expectedSpecInfo, 'AbsTol', 1e-10);
        end

         % function testCalculateAtom(testCase)
         %     
         % end

      
    end

end