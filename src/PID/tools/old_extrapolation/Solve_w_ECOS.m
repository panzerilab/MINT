classdef Solve_w_ECOS < handle
%%% BROJA_2PID.py -- Python module
%%%
%%% BROJA_2PID: Bertschinger-Rauh-Olbrich-Jost-Ay (BROJA) bivariate Partial Information Decomposition
%%% https://github.com/Abzinger/BROJA_2PID
%%% (c) Abdullah Makkeh, Dirk Oliver Theis
%%% Permission to use and modify with proper attribution
%%% (Apache License version 2.0)
%%%
%%% Information about the algorithm, documentation, and examples are here:
%%% @Article{makkeh-theis-vicente:pidOpt:2017,
%%%          author =       {Makkeh, Abdullah and Theis, Dirk Oliver and Vicente, Raul},
%%%          title =        {BROJA-2PID: A cone programming based Partial Information Decomposition estimator},
%%%          journal =      {jo},
%%%          year =         2017,
%%%          key =       {key},
%%%          volume =    {vol},
%%%          number =    {nr},
%%%          pages =     {1--2}
%%% }
%%% Please cite this paper when you use this software (cf. README.md)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%
%%% ECOS's exp cone: (r,p,q)   w/   q>0  &  exp(r/q) <= p/q
%%% Translation:     (0,1,2)   w/   2>0  &  0/2      <= ln(1/2)
%%%
%%% Methods:
%%%   - Solve_w_ECOS: Constructor method
%%%   - create_model: Creates the optimization model
%%%   - r_vidx: Computes the index for variable r_i
%%%   - p_vidx: Computes the index for variable p_i
%%%   - q_vidx: Computes the index for variable q_i
%%%   - solve: Solves the optimization problem
%%%   - entropy_X: Computes the entropy of variable X
%%%   - entropy_X2: Computes the entropy of variable X using precomputed table
%%%   - condenttropy: Computes conditional entropy
%%%   - addqdistr: Adds quantum distribution to the triplets
%%%   - condentropy2: Computes conditional entropy using precomputed table
%%%   - condentropy_orig: Computes original conditional entropy
%%%   - create_qxyz: Creates q_xyz table
%%%   - condYmutinf: Computes conditional mutual information for Y
%%%   - condZmutinf: Computes conditional mutual information for Z
%%%   - dual_value: Computes the dual value
%%%   - check_feaseability: Checks the feasibility of the solution
    
    properties
        ecos_kwargs         %Structure for ECOS parameters
        verbose             %Controls the output of debug information
        c                   %Parameters for the convex optimization problem
        G
        h
        dims                %Structure for the dimensions of the optimization problem
        A                   %Parameters for the constraints of the optimization problem
        b
        sol_rpq             %Attributes for storing solution information
        sol_slack
        sol_lambda
        sol_mu
        sol_info
        solution            %Structure for the solution of the optimization problem
        b_xy                %Probability density functions for the conditions of PID
        b_xz
        NX                  %Number of states for variables X, Y, Z
        NY
        NZ
        triplet_nonzero     %Table of non-empty triplets (X, Y, Z) and their associated indices
        n                   %Number of triplets and constraints
        m
        table_pxy           % Tables for the marginal distributions
        table_pxz
        table_pyz
        table_px
        Eqn_list            % List of equations
    end

    methods
        function obj = Solve_w_ECOS(marg_xy, marg_xz)
            % (c) Abdullah Makkeh, Dirk Oliver Theis
            % Permission to use and modify under Apache License version 2.0

            % ECOS parameters
            obj.ecos_kwargs = struct();
            obj.verbose = false;

            % Data for ECOS
            % obj.c = [];
            % obj.G = [];
            % obj.h = [];
            % obj.dims = struct();
            % obj.A = [];
            % obj.b = [];


            % ECOS result
            obj.sol_rpq = [];
            obj.sol_slack = [];
            obj.sol_lambda = [];
            obj.sol_mu = [];
            obj.sol_info = [];
            obj.solution= struct;

            % Probability density function data
            % marg_xy = sum(prob_xyz, 3);
            % marg_xz = sum(prob_xyz, 2);
            % marg_xy = reshape(marg_xy, [size(marg_xy, 1) size(marg_xy, 2)]);
            % marg_xz = reshape(marg_xz, [size(marg_xz, 1) size(marg_xz, 3)]);
            obj.b_xy = marg_xy;
            obj.b_xz = marg_xz;
            obj.NX = size(obj.b_xy, 1);
            obj.NY = size(obj.b_xy, 2);
            obj.NZ = size(obj.b_xz, 2);


            % disp('NX')
            % disp(obj.NX)
            % disp('NY')
            % disp(obj.NY)
            % disp('NZ')
            % disp(obj.NZ)

            %n = number of  (x, y, z) triplets with non-zero probability
            %m = number of (x, y) pairs with non-zero probability+number of (x, z) pairs with non-zero probability
            assert(ndims(obj.b_xy) == 2)
            assert(ndims(obj.b_xz) == 2)
            [nzXy, nzY, pxy] = find(obj.b_xy);
            [nzXz, nzZ, pxz] = find(obj.b_xz);
            if isrow(nzXy)
                nzXy = nzXy';
                nzY = nzY';
                pxy = pxy';
            end

            if isrow(nzXz)
                nzXz = nzXz';
                nzZ = nzZ';
                pxz = pxz';
            end

            try
                obj.table_pxy = array2table([nzXy nzY pxy], 'VariableNames',{'X', 'Y', 'pxy'} );
                obj.table_pxz = array2table([nzXz nzZ pxz], 'VariableNames',{'X', 'Z', 'pxz'} );
            catch
                disp('nzXy')
                disp(nzXy)
                disp('nzY')
                disp(nzY)
                disp('pxy')
                disp(pxy)
                disp('nzXz')
                disp(nzXz)
                disp('nzZ')
                disp(nzZ)
                disp('pxz')
                disp(pxz)
                disp('obj.table_pxy')
                disp(obj.table_pxy)
                disp('obj.table_pxz')
                disp(obj.table_pxz)
            end
            
            try
                obj.triplet_nonzero = innerjoin(obj.table_pxy, obj.table_pxz);
                % obj.triplet_nonzero = obj.triplet_nonzero(:, {'X', 'Y', 'Z'});
            catch
                disp('obj.table_pxy')
                disp(obj.table_pxy)
                disp('obj.table_pxz')
                disp(obj.table_pxz)
            end
            % gx  = groupsummary(obj.table_pxy,{'X'}, 'sum');
            % gx  = gx(:,{'X','sum_pxy'});
            % gx  = renamevars(gx,'sum_pxy',  "p_x");
            % gy  = groupsummary(obj.table_pxy,{'Y'}, 'sum');
            % gy  = gy(:,{'Y','sum_pxy'});
            % gy  = renamevars(gy,'sum_pxy',  "p_y");
            % 
            % obj.triplet_nonzero = innerjoin(obj.triplet_nonzero, gx);
            % obj.triplet_nonzero = innerjoin(obj.triplet_nonzero, gy);
            
            obj.triplet_nonzero.index = (1:height(obj.triplet_nonzero)).';
            obj.n = height(obj.triplet_nonzero); % size(obj.triplet_nonzero,1);
            obj.m = height(obj.table_pxy) + height(obj.table_pxz);% nnz(obj.b_xy) + nnz(obj.b_xz);
            % obj.table_px = groupsummary(obj.table_pxy,{'X'},'sum');
            % obj.table_px = renamevars(obj.table_px, "sum_pxy",  "p_x");
            % marg_yz = sum(prob_xyz, 1);
            % marg_yz = reshape(marg_yz, [size(marg_yz, 2) size(marg_yz, 3)]);
            % [nzYz, nzZ, pyz] = find(marg_yz);
            % obj.table_pyz = array2table([nzYz nzZ pyz], 'VariableNames',{'Y', 'Z', 'pyz'} );

        end

        function create_model(obj)
            % (c) Abdullah Makkeh, Dirk Oliver Theis
            % Permission to use and modify under Apache License version 2.0

            n_vars = 3 * obj.n;
            n_cons = obj.n + obj.m;

            % Create the equations: Ax = b
            obj.b = zeros(n_cons, 1,'double');

            % Eqn = zeros(2 * obj.n + obj.m, 1,'double');
            % Var = zeros(obj.n, 1,'double');
            % Coeff = zeros(n_cons, 1,'double');

            Eqn = [];
            Var = [];
            Coeff = [];

            % The q-p coupling equations: q_{*yz} - p_{xyz} = 0
            % disp('The q-p coupling equations: q_{*yz} - p_{xyz} = 0')
            % disp('n')
            % disp(obj.n)
            % disp('m')
            % disp(obj.m)
            eqn = 0;
            for i = 1:obj.n
                xyz = obj.triplet_nonzero(i,:);
                eqn = eqn+1;
                p_var = obj.p_vidx(i);
                Eqn = [Eqn; eqn];
                % Eqn(i) = i;
                Var = [Var; p_var];
                % Var(i) = p_var;
                Coeff = [Coeff; -1];
                % Coeff(i) = -1;

                %x = xyz{1, 1};
                y = xyz{1, 'Y'}; %xyz{1, 2};
                z = xyz{1, 'Z'}; %xyz{1, 3};

                select_idx = obj.triplet_nonzero.Y == y & obj.triplet_nonzero.Z == z;
                select_tab = obj.triplet_nonzero(select_idx,:);
                if height(select_tab)>0
                    select_tab.qid = obj.q_vidx(select_tab.index);
                    Eqn = [Eqn; repmat(eqn, height(select_tab), 1)];
                    Var = [Var; select_tab.qid];
                    Coeff = [Coeff; ones(height(select_tab), 1)];
                end

            end

            % running number
            %eqn = obj.n;
            % disp('last index of the coupling pq equations')
            % disp(eqn)
            % The xy marginals q_{xy*} = b^y_{xy}
            % disp('The xy marginals q_{xy*} = b^y_{xy}')
            for row=1:height(obj.table_pxy)
                eqn = eqn + 1;
                selec_triplets = innerjoin(obj.triplet_nonzero, obj.table_pxy(row,:));
                % selec_triplets = splitapply(@sum,Weight,G)
                if height(selec_triplets)>0
                    selec_triplets.qid = obj.q_vidx(selec_triplets.index);
                    Eqn = [Eqn; repmat(eqn, height(selec_triplets), 1)];
                    Var = [Var; selec_triplets.qid];
                    Coeff = [Coeff; ones(height(selec_triplets), 1)];
                end

                obj.b(eqn) = obj.b_xy(obj.table_pxy{row,1}, obj.table_pxy{row,2});
            end

            % disp('last index of the xy marginals')
            % disp(eqn)

            % The xz marginals q_{xz*} = b^z_{xz}
            % disp('The xz marginals q_{xz*} = b^z_{xz}')
            for row=1:height(obj.table_pxz)
                eqn = eqn + 1;
                selec_triplets = innerjoin(obj.triplet_nonzero, obj.table_pxz(row,:));
                if height(selec_triplets)>0
                    selec_triplets.qid = obj.q_vidx(selec_triplets.index);
                    Eqn = [Eqn; repmat(eqn, height(selec_triplets), 1)];
                    Var = [Var; selec_triplets.qid];
                    Coeff = [Coeff; ones(height(selec_triplets), 1)];
                end

                obj.b(eqn) = obj.b_xz(obj.table_pxz{row,1}, obj.table_pxz{row,2});
            end
            % disp('last index of the xz marginals')
            % disp(eqn)
            obj.Eqn_list = Eqn;
            % disp('last index of the Eqn')
            % disp(Eqn(end))
            obj.A = sparse(Eqn, Var, Coeff, n_cons, n_vars);
            % disp('size(obj.A)')
            % disp(size(obj.A))
            % Generalized ieqs: gen.nneg of the variable triples (r_i,q_i,p_i), i=0,dots,n-1:
            Ieq = [];
            IeqVar = [];
            Coeff = [];
            for i = 1:obj.n
                %xyz = obj.trip_of_idx{i};
                r_var = obj.r_vidx(i);
                q_var = obj.q_vidx(i);
                p_var = obj.p_vidx(i);

                Ieq = [Ieq; length(Ieq) + 1; length(Ieq) + 2; length(Ieq) + 3];
                IeqVar = [IeqVar; r_var; p_var; q_var];
                Coeff = [Coeff; -1; -1; -1];
            end

            obj.G = sparse(Ieq, IeqVar, Coeff, n_vars, n_vars);
            obj.h = zeros(n_vars, 1);
            obj.dims.e = obj.n;

            % Objective function:
            obj.c = zeros(n_vars, 1);
            for i = 1:obj.n
                obj.c(obj.r_vidx(i)) = -1;
            end
        end

        function r = r_vidx(obj, i)
            r = 3 * (i-1) + 1;
        end

        function p = p_vidx(obj, i)
            p = 3 * (i-1) + 2;
        end

        function q = q_vidx(obj, i)
            q = 3 * (i-1) + 3;
        end

        function result = solve(obj)
            % (c) Abdullah Makkeh, Dirk Oliver Theis
            % Permission to use and modify under Apache License version 2.0
            % for cond[]mutinf computation below

            if ~isempty(obj.verbose)
                obj.ecos_kwargs.verbose = obj.verbose;
            end
            obj.ecos_kwargs.feastol = 1E-7;
            obj.ecos_kwargs.abstol = 1E-6;
            obj.ecos_kwargs.reltol = 1E-6;

            % obj.solution =
            [x,y,info,s,z] = ecos(obj.c, obj.G, obj.h, obj.dims, obj.A, obj.b, obj.ecos_kwargs);

            if size(x,1)> 0
                obj.sol_rpq = x;
                obj.sol_slack = s;
                obj.sol_lambda = y;
                obj.sol_mu = z;
                obj.sol_info = info;
                result = 'success';
            else
                % "x" not in struct solution
                result = 'x not in struct solution -- No Solution Found!!!';
            end
        end

        function result = entropy_X(obj, pdf)
            [indX, indYZ, pvalues] = find(pdf);
            [indY, indZ] = ind2sub([obj.NY obj.NZ], indYZ);

            table_pxyz = array2table([indX, indY, indZ, pvalues],'VariableNames',{'X','Y','Z', 'p'} );

            g = groupsummary(table_pxyz,{'X'},'sum');
            g.('entr') = - g.('sum_p').*log(g.('sum_p'));
            result = sum(g{:,'entr'},1);
        end

        function result = entropy_X2(obj)
            g = groupsummary(obj.table_pxy,{'X'},'sum');
            g.('entr') = - g.('sum_pxy').*log(g.('sum_pxy'));
            result = sum(g{:,'entr'},1);
            gx  = g(:,{'X','sum_pxy'});
            gx.Properties.VariableNames{'sum_pxy'} = 'p_x';
            obj.triplet_nonzero = innerjoin(obj.triplet_nonzero, gx);
        end

        function result=condenttropy(obj)
            result = 0;
            for yv =1:obj.NY
                for zv = 1:obj.NZ
                    marg_x = 0;
                    indexes = find(obj.triplet_nonzero.Y == yv & obj.triplet_nonzero.Z == zv);s
                    if size(indexes,1) ~= 0
                        for i = 1:size(indexes,1)
                            q_id = obj.q_vidx(indexes(i));
                            marg_x = marg_x + max(0,self.sol_rpq(q_id));
                        end
                        for i = 1:size(indexes,1)
                            q_id = obj.q_vidx(indexes(i));
                            q = self.sol_rpq(q_id);
                            if q > 0
                                result = result - q*log(q/marg_x);
                            end
                        end
                    end
                end
            end
        end

        function addqdistr(obj)
            obj.triplet_nonzero.qid = 3 * obj.triplet_nonzero.index;
            obj.triplet_nonzero.q = obj.sol_rpq(obj.triplet_nonzero.qid);
        end

        function create_qxyz(obj)
            obj.triplet_nonzero.qpos = obj.triplet_nonzero.q;
            obj.triplet_nonzero.qpos(obj.triplet_nonzero.qpos<0) = 0;
            
            gyz = groupsummary(obj.triplet_nonzero,{'Y', 'Z'}, 'sum');
            gyz = gyz(:,{'Y', 'Z','sum_qpos'});
            gyz.Properties.VariableNames{'sum_qpos'} = 'q_yz';

            gy  = groupsummary(obj.triplet_nonzero,{'Y'}, 'sum');
            gy  = gy(:,{'Y','sum_qpos'});
            gy.Properties.VariableNames{'sum_qpos'} = 'q_y';
           
            gz  = groupsummary(obj.triplet_nonzero,{'Z'}, 'sum');
            gz  = gz(:,{'Z','sum_qpos'});
            gz.Properties.VariableNames{'sum_qpos'} = 'q_z';

            obj.triplet_nonzero = innerjoin(obj.triplet_nonzero, gyz);
            obj.triplet_nonzero = innerjoin(obj.triplet_nonzero, gy);
            obj.triplet_nonzero = innerjoin(obj.triplet_nonzero, gz);
        end

        function result=condentropy2(obj)
            % gyz = groupsummary(obj.triplet_nonzero,{'Y', 'Z'}, 'sum');
            % gyz = gyz(:,{'Y', 'Z','sum_qpos'});
            % gyz = renamevars(gyz,"sum_qpos",  "qpos_yz");
            % obj.triplet_nonzero = innerjoin(obj.triplet_nonzero, gyz);
            mask = obj.triplet_nonzero.qpos > 0;
            condentr = -obj.triplet_nonzero.qpos(mask) .* log(obj.triplet_nonzero.qpos(mask) ./ obj.triplet_nonzero.q_yz(mask));
            result = sum(condentr);
        end

        function result = condentropy_orig(obj, pdf)
            [indX, indYZ, pvalues] = find(pdf);
            [indY, indZ] = ind2sub([obj.NY obj.NZ], indYZ);

            if isrow(indX)
                indX = indX';
                indY = indY';
                indZ = indZ';
                pvalues = pvalues';
            end

            table_pxyz = array2table([indX, indY, indZ, pvalues],'VariableNames',{'X','Y','Z', 'p'} );

            g = groupsummary(table_pxyz,{'Y', 'Z'},'sum');
            tmp_table = innerjoin(table_pxyz, g);
            tmp_table.('condentrorig') = -tmp_table.('p') .* log(tmp_table.('p') ./ tmp_table.('sum_p'));
            % condentrorig = -obj.triplet_nonzero.p .* log(obj.triplet_nonzero.p ./ obj.triplet_nonzero.q_yz);
            result = sum(tmp_table{:,'condentrorig'},1);%sum(condentrorig);
        end
        
        function result = condentropy_orig2(obj)
            tmp_table.('condentrorig') = -obj.triplet_nonzero.('p') .* log(obj.triplet_nonzero.('p') ./ obj.triplet_nonzero.('sum_p'));
            % condentrorig = -obj.triplet_nonzero.p .* log(obj.triplet_nonzero.p ./ obj.triplet_nonzero.q_yz);
            result = sum(tmp_table{:,'condentrorig'},1);%sum(condentrorig);
        end
        
        function result = condYmutinf(obj)
            tmp_table = innerjoin(obj.triplet_nonzero, obj.table_pxy);
            % obj.triplet_nonzero = sortrows(obj.triplet_nonzero,"index","ascend");
            tmp_table = sortrows(tmp_table,"index","ascend");
            % obj.triplet_nonzero.condYmutinf = obj.triplet_nonzero.q .* log(obj.triplet_nonzero.q .* obj.triplet_nonzero.q_y ./ (tmp_table.pxy .* obj.triplet_nonzero.q_yz));
            mask = tmp_table.qpos > 0;
            condYmutinf = tmp_table.qpos(mask) .* log(tmp_table.qpos(mask) .* tmp_table.q_y(mask) ./ (tmp_table.pxy(mask) .* tmp_table.q_yz(mask)));
            result = sum(condYmutinf, 'all'); %sum(tmp_table{:,'condYmutinf'},1);
        end

        function result = condYmutinf2(obj)
            % obj.triplet_nonzero = sortrows(obj.triplet_nonzero,"index","ascend");
            % obj.triplet_nonzero.condYmutinf = obj.triplet_nonzero.q .* log(obj.triplet_nonzero.q .* obj.triplet_nonzero.q_y ./ (tmp_table.pxy .* obj.triplet_nonzero.q_yz));
            mask = obj.triplet_nonzero.qpos > 0;
            condYmutinf = obj.triplet_nonzero.qpos(mask) .* log(obj.triplet_nonzero.qpos(mask) .* obj.triplet_nonzero.q_y(mask) ./ (obj.triplet_nonzero.pxy(mask) .* obj.triplet_nonzero.q_yz(mask)));
            result = sum(condYmutinf, 'all'); %sum(tmp_table{:,'condYmutinf'},1);
        end
        function result = condZmutinf(obj)
            tmp_table = innerjoin(obj.triplet_nonzero, obj.table_pxz);
            % obj.triplet_nonzero = sortrows(obj.triplet_nonzero,"index","ascend");
            tmp_table = sortrows(tmp_table,"index","ascend");
            % obj.triplet_nonzero.condZmutinf = obj.triplet_nonzero.q .* log(obj.triplet_nonzero.q .* obj.triplet_nonzero.q_z ./ (tmp_table.pxz .* obj.triplet_nonzero.q_yz));
            mask = tmp_table.qpos > 0;
            condZmutinf = tmp_table.qpos(mask) .* log(tmp_table.qpos(mask) .* tmp_table.q_z(mask) ./ (tmp_table.pxz(mask) .* tmp_table.q_yz(mask)));
            result = sum(condZmutinf, 'all'); % sum(tmp_table{:,'condZmutinf'},1);
        end

        function result = condZmutinf2(obj)
            % obj.triplet_nonzero = sortrows(obj.triplet_nonzero,"index","ascend");
            % obj.triplet_nonzero.condZmutinf = obj.triplet_nonzero.q .* log(obj.triplet_nonzero.q .* obj.triplet_nonzero.q_z ./ (tmp_table.pxz .* obj.triplet_nonzero.q_yz));
            mask = obj.triplet_nonzero.qpos > 0;
            condZmutinf = obj.triplet_nonzero.qpos(mask) .* log(obj.triplet_nonzero.qpos(mask) .* obj.triplet_nonzero.q_z(mask) ./ (obj.triplet_nonzero.pxz(mask) .* obj.triplet_nonzero.q_yz(mask)));
            result = sum(condZmutinf, 'all'); % sum(tmp_table{:,'condZmutinf'},1);
        end

        function result = dual_value(obj)
            result = -dot(obj.sol_lambda, obj.b);
        end

        function [primal_infeasibility, dual_infeasibility] = check_feasibility(obj)
            max_q_negativity = max(0, -min(obj.triplet_nonzero.q));
            obj.triplet_nonzero.max_q_positive = max(0, obj.triplet_nonzero.q);
            max_violation_of_eqn = 0;
            gxy = groupsummary(obj.triplet_nonzero,{'X', 'Y'}, 'sum');
            gxy = gxy(:,{'X', 'Y','sum_max_q_positive'});
            gxy.Properties.VariableNames{'sum_max_q_positive'} = 'sum_max_qxy_positive';
            tmp_table = innerjoin(obj.table_pxy,gxy);
            tmp_table.eqnviolation = abs(tmp_table.pxy - tmp_table.sum_max_qxy_positive);
            max_violation_of_eqn = max(max_violation_of_eqn, max(tmp_table.eqnviolation));

            gxz = groupsummary(obj.triplet_nonzero,{'X', 'Z'}, 'sum');
            gxz = gxz(:,{'X', 'Z','sum_max_q_positive'});
            gxz.Properties.VariableNames{'sum_max_q_positive'} = 'sum_max_qxz_positive';
            tmp_table = innerjoin(obj.table_pxz,gxz);
            tmp_table.eqnviolation = abs(tmp_table.pxz - tmp_table.sum_max_qxz_positive);
            max_violation_of_eqn = max(max_violation_of_eqn, max(tmp_table.eqnviolation));
            primal_infeasibility = max(max_violation_of_eqn, max_q_negativity);
            dual_infeasibility = 0;
            % obj.triplet_nonzero.mu_yz = -self.sol_lambda[xy_idx] - self.sol_lambda[xz_idx] - mu_yz[(y,z)] -ln(-self.sol_lambda[i]) - 1

        end
    end
end

