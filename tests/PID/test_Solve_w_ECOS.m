classdef test_Solve_w_ECOS < matlab.unittest.TestCase

    methods (Test)
        function testConstructor(testCase)
            rng("default")
            marg_xy = rand(100,100);
            marg_xy = marg_xy / max(marg_xy, [], 'all');
            marg_xz = rand(100,100);
            marg_xz = marg_xz / max(marg_xz, [], 'all');
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            testCase.verifyClass(myObject, 'Solve_w_ECOS');
        end

        function calculate_rand_value(testCase)
            rng("default")
            marg_xy = sprand(30,10,0.2); 
            marg_xy = marg_xy / sum(marg_xy, 'all');
            marg_xz = sprand(30,20,0.2); 
            marg_xz = marg_xz / sum(marg_xz, 'all');
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
        end

        function calculate_and_gate(testCase)
            marg_xy = zeros(2);
            marg_xy(1,1) = .5;
            marg_xy(1,2) = .25;
            marg_xy(2,2) = .25;
            marg_xz = zeros(2);
            marg_xz(1,1) = .5;
            marg_xz(1,2) = .25;
            marg_xz(2,2) = .25;

            p_xyz = zeros(2, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end

            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
        end

        function calculate_copy_gate(testCase)
            p_xyz = zeros(2, 2, 2);
            p_xyz(1,1,1) = .25;
            p_xyz(2,2,1) = .25;
            p_xyz(1,1,2) = .25;
            p_xyz(2,2,2) = .25;

            marg_xy = squeeze(sum(p_xyz, 2));
            marg_xz = squeeze(sum(p_xyz, 3));

            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            assert(isequal(result, 'success'));
        end

        function calculate_entropy_X(testCase)
            rng("default")
            prob_xyz = rand(10, 15, 5);
            prob_xyz = prob_xyz / sum(prob_xyz, 'all');
            prob_xyz = permute(prob_xyz, [3 2 1]);
            marg_xy = squeeze(sum(prob_xyz, 2));
            marg_xz = squeeze(sum(prob_xyz, 3));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            entropy_X2 = myObject.entropy_X2();
            assert(isequal(result, 'success'));
            assert(abs(entropy_X2 - 1.6093) < 0.001);
        end

        function calculate_entropy_X_and_gate(testCase)
            marg_xy = zeros(2);
            marg_xy(1,1) = .5;
            marg_xy(1,2) = .25;
            marg_xy(2,2) = .25;
            marg_xz = zeros(2);
            marg_xz(1,1) = .5;
            marg_xz(1,2) = .25;
            marg_xz(2,2) = .25;
            p_xyz = zeros(2, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end
            p_xyz = permute(p_xyz, [3 2 1]);
            marg_xy = squeeze(sum(p_xyz, 2));
            marg_xz = squeeze(sum(p_xyz, 3));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();            
            entropy_X2 = myObject.entropy_X2();
            assert(isequal(result, 'success'));
            assert(abs(entropy_X2 - 0.6931) < 0.001);
        end

        function calculate_addqdistr(testCase)
            rng("default")
            marg_xy = zeros(2);
            marg_xy(1,1) = .5;
            marg_xy(1,2) = .25;
            marg_xy(2,2) = .25;
            marg_xz = zeros(2);
            marg_xz(1,1) = .5;
            marg_xz(1,2) = .25;
            marg_xz(2,2) = .25;
            p_xyz = zeros(2, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end
            p_xyz = permute(p_xyz, [3 2 1]);
            marg_xy = squeeze(sum(p_xyz, 2));
            marg_xz = squeeze(sum(p_xyz, 3));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            assert(isequal(result, 'success'));
            myObject.addqdistr()
            assert(~isempty(myObject.triplet_nonzero.qid))
        end

        function calculate_large_distribution(testCase)
            rng("default")
            n = 4;
            p_xyz = ones(n, n, n);
            p_xyz = p_xyz / sum(p_xyz, 'all');
            marg_xy = squeeze(sum(p_xyz, 2));
            marg_xz = squeeze(sum(p_xyz, 3));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X2();
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf();
            condZmutinf   = myObject.condZmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf* bits;
            uiz = condZmutinf* bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0,0,0,0];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 1.3863) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end

        function calculate_condentropy_orig(testCase)
            rng("default")
            prob_xyz = rand(20,15, 5);
            prob_xyz = prob_xyz / sum(prob_xyz, 'all');
            marg_xy = squeeze(sum(prob_xyz, 2));
            marg_xz = squeeze(sum(prob_xyz, 3));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model();
            result = myObject.solve();
            myObject.addqdistr();
            condent_orig = myObject.condentropy_orig(prob_xyz);
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(condent_orig - 2.8204) < 0.001);
        end

        function calculate_condYmutinf(testCase)
            rng("default")
            prob_xyz = rand(5,10, 5); 
            prob_xyz = prob_xyz / sum(prob_xyz, 'all');
            marg_xy = squeeze(sum(prob_xyz, 2));
            marg_xz = squeeze(sum(prob_xyz, 3));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr();
            myObject.create_qxyz();
            conYmutinf = myObject.condYmutinf();
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.q_yz));
            assert(~isempty(myObject.triplet_nonzero.q_y));
            assert(~isempty(myObject.triplet_nonzero.q_z));
            assert(~isempty(myObject.triplet_nonzero.qid));
            assert(abs(conYmutinf - 0.0081) < 0.001);
        end

        function calculate_condZmutinf(testCase)
            rng("default")
            prob_xyz = rand(20,15, 5); 
            prob_xyz = prob_xyz / sum(prob_xyz, 'all');
            marg_xy = squeeze(sum(prob_xyz, 2));
            marg_xz = squeeze(sum(prob_xyz, 3));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model();
            result = myObject.solve();
            myObject.addqdistr();
            myObject.create_qxyz();
            condZmutinf = myObject.condZmutinf();
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.q_yz));
            assert(~isempty(myObject.triplet_nonzero.q_y));
            assert(~isempty(myObject.triplet_nonzero.q_z));
            assert(~isempty(myObject.triplet_nonzero.qid));
            assert(abs(condZmutinf - 0.0064) < 0.001);
        end

        function calculate_CI_copy_gate(testCase)        
            marg_xy = zeros(2);
            marg_xy(1,1) = .5;
            marg_xy(2,2) = .5;
            marg_xz = ones(2)/4;
            p_xyz = zeros(2, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model();
            result = myObject.solve();
            myObject.addqdistr();
            myObject.create_qxyz();
            entropy_X     = myObject.entropy_X(p_xyz);
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condZmutinf();
            condZmutinf   = myObject.condYmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf* bits;
            uiz = condZmutinf* bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0,1,0,0];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 0.6931) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end

        
        function calculate_CI_and_gate(testCase)
            marg_xy = zeros(2);
            marg_xy(1,1) = .5;
            marg_xy(1,2) = .25;
            marg_xy(2,2) = .25;
            marg_xz = zeros(2);
            marg_xz(1,1) = .5;
            marg_xz(1,2) = .25;
            marg_xz(2,2) = .25;
            p_xyz = zeros(2, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X2();
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf();
            condZmutinf   = myObject.condZmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf * bits;
            uiz = condZmutinf * bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0.311,0,0,0.229];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 0.5623) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end

        function calculate_CI_or_gate(testCase)
            marg_xy = zeros(2);
            marg_xy(1,1) = .25;
            marg_xy(2,1) = .25;
            marg_xy(2,2) = .5;
            marg_xz = marg_xy;
            p_xyz = zeros(2, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X(p_xyz);
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf();
            condZmutinf   = myObject.condZmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf* bits;
            uiz = condZmutinf* bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0.311,0,0,0.229];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 0.5623) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end

        function calculate_CI_reducedor_gate(testCase)
            p_xyz = zeros(2, 2, 2);
            p_xyz(1,1,1)=.5;
            p_xyz(1,2,2)=.25;
            p_xyz(2,1,2)=.25;
            marg_xy = sum(p_xyz, 3);
            marg_xz = squeeze(sum(p_xyz, 2));
            p_x = sum(marg_xz,2);
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X(p_xyz);
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf2();
            condZmutinf   = myObject.condZmutinf2();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf* bits;
            uiz = condZmutinf* bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0.122,0.188,0,0.5];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 0.5623) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.01);
        end


        function calculate_CI_xor_gate(testCase)
            marg_xy = ones(2)/4;
            marg_xz = marg_xy;
            p_xyz = zeros(2, 2, 2);

            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            p_xyz(:,:,1) = [.25 0; 0 .25];
            p_xyz(:,:,2) = [0 .25; .25 0];

            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X2();
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf();
            condZmutinf   = myObject.condZmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf* bits;
            uiz = condZmutinf* bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0,0,0,1];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 0.6931) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end
     
        function calculate_CI_copyboth_gate(testCase)
            marg_xy = zeros(4,2);
            marg_xy(1,1) = .25;
            marg_xy(2,1) = .25;
            marg_xy(3,2) = .25;
            marg_xy(4,2) = .25;
            marg_xz = zeros(4,2);
            marg_xz(1,1) = .25;
            marg_xz(2,2) = .25;
            marg_xz(3,1) = .25;
            marg_xz(4,2) = .25;

            p_xyz = zeros(4, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X2();
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf();
            condZmutinf   = myObject.condZmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf* bits;
            uiz = condZmutinf* bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0,1,1,0];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 1.3863) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end
       
        function calculate_CI_copyboth_gate4(testCase)
            marg_xy = zeros(4,2);
            marg_xy(1,1) = .25;
            marg_xy(2,1) = .25;
            marg_xy(3,2) = .25;
            marg_xy(4,2) = .25;
            marg_xz = zeros(4,2);
            marg_xz(1,1) = .25;
            marg_xz(2,2) = .25;
            marg_xz(3,1) = .25;
            marg_xz(4,2) = .25;

            p_xyz = zeros(4, 2, 2);
            p_x = sum(marg_xz,2);
            % Assuming conditional independence: p(y, z | x) = p(y | x) * p(z | x)
            for ix = 1:2
                for iy=1:2
                    for iz = 1:2
                        p_xyz(ix, iy, iz) = marg_xy(ix, iy)*marg_xz(ix, iz)/p_x(ix);
                    end
                end
            end
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X2();
            condent       = myObject.condentropy2();
            condent_orig  = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf();
            condZmutinf   = myObject.condZmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf* bits;
            uiz = condZmutinf* bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [0,1,1,0];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 1.3863) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end

        function calculate_CI_xorand_gate(testCase)
            % Test the that we reach to a solution
            p_xyz = zeros(4, 2, 2);
            p_xyz(1,1,1) = .25;
            p_xyz(3,2,1) = .25;
            p_xyz(2,2,2) = .25;
            p_xyz(3,1,2) = .25;

            marg_xy = sum(p_xyz,3);
            marg_xz = squeeze(sum(p_xyz,2));

            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            bits = 1/log(2);
            entropy_X     = myObject.entropy_X2()* bits;
            condent       = myObject.condentropy2()* bits;
            condent_orig = myObject.condentropy_orig(p_xyz)* bits;
            condYmutinf   = myObject.condYmutinf()* bits;
            condZmutinf   = myObject.condZmutinf()* bits;
            
            si  = (entropy_X  - condent - condZmutinf - condYmutinf);
            uiy = condYmutinf;
            uiz = condZmutinf;
            ci  = (condent - condent_orig);
            pid_v = [si uiy uiz ci];
            expected_pid_v = [.5,0.,0.,1.];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end

        function calculate_CI_redundant_gate(testCase)
            p_xyz = zeros(2, 2, 2);
            for ix = 1:2
                p_xyz(ix,ix,ix)=.5;
            end
            marg_xy = sum(p_xyz,3);
            marg_xz = squeeze(sum(p_xyz,2));
            myObject = Solve_w_ECOS(marg_xy, marg_xz);
            myObject.create_model()
            result = myObject.solve();
            myObject.addqdistr()
            myObject.create_qxyz()
            entropy_X     = myObject.entropy_X2();
            condent       = myObject.condentropy2();
            condent_orig = myObject.condentropy_orig(p_xyz);
            condYmutinf   = myObject.condYmutinf();
            condZmutinf   = myObject.condZmutinf();
            bits = 1/log(2);
            si  = (entropy_X  - condent - condZmutinf - condYmutinf) * bits;
            uiy = condYmutinf * bits;
            uiz = condZmutinf * bits;
            ci  = (condent - condent_orig)* bits;
            pid_v = [si uiy uiz ci];
            expected_pid_v = [1,0,0,0];
            assert(isequal(result, 'success'));
            assert(~isempty(myObject.triplet_nonzero.qid))
            assert(abs(entropy_X - 0.6931) < 0.001);
            assert(abs(sum(pid_v - expected_pid_v)) < 0.001);
        end
   
    end
end
