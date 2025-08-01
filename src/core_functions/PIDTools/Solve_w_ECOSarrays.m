classdef Solve_w_ECOSarrays < handle
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
        ecos_kwargs         % Structure for ECOS parameters
        verbose             % Controls the output of debug information
        c                   % Parameters for the convex optimization problem
        G
        h
        dims                % Structure for the dimensions of the optimization problem
        A                   % Parameters for the constraints of the optimization problem
        b
        sol_rpq             % Attributes for storing solution information
        sol_slack
        sol_lambda
        sol_mu
        sol_info
        solution            % Structure for the solution of the optimization problem
        b_xy                % Probability density functions for the conditions of PID
        b_xz
        NX                  % Number of states for variables X, Y, Z
        NY
        NZ
        triplet_nonzero     % Matrix of non-empty triplets (X, Y, Z) and their associated indices
        nzXYZx
        nzXYZyz
        n                   % Number of triplets and constraints
        m
        pxy                 % Matrices for the marginal distributions
        pxz
        nzXy
        nzXz
        nzY
        nzZ
        qyz
        qy
        qz
        triplet_nonzeroqid
        triplet_nonzeroq
        triplet_nonzeroqpos
        Eqn_list            % List of equations
        Var_list
        Coeff_list
    end

    methods
        function obj = Solve_w_ECOSarrays(marg_xy, marg_xz)
            % (c) Abdullah Makkeh, Dirk Oliver Theis
            % Permission to use and modify under Apache License version 2.0

            % ECOS parameters
            obj.ecos_kwargs = struct();
            obj.verbose = false;

            % Data for ECOS
            obj.c = [];
            obj.G = [];
            obj.h = [];
            obj.dims = struct();
            obj.A = [];
            obj.b = [];
            obj.Eqn_list = [];
            obj.Var_list = [];
            obj.Coeff_list = [];

            % ECOS result
            obj.sol_rpq = [];
            obj.sol_slack = [];
            obj.sol_lambda = [];
            obj.sol_mu = [];
            obj.sol_info = [];
            obj.solution = struct;

            % Probability density function data
            obj.b_xy = marg_xy;
            obj.b_xz = marg_xz;
            obj.NX = size(obj.b_xy, 1);
            obj.NY = size(obj.b_xy, 2);
            obj.NZ = size(obj.b_xz, 2);

            [obj.nzXy, obj.nzY, obj.pxy] = find(obj.b_xy);
            [obj.nzXz, obj.nzZ, obj.pxz] = find(obj.b_xz);

            % digits(50)
            % [nonzeros_per_rowY, ~] =  histcounts(obj.nzXy,'BinMethod','integers'); %hist(obj.nzXy,unique(obj.nzXy));
            % [nonzeros_per_rowZ, ~] =  histcounts(obj.nzXz,'BinMethod','integers'); %hist(obj.nzXz,unique(obj.nzXz));

            nonzeros_per_rowY = sum(obj.b_xy>0,2);
            nonzeros_per_rowZ = sum(obj.b_xz>0,2);
            obj.n =  dot(nonzeros_per_rowY,nonzeros_per_rowZ);
            obj.m = size(obj.pxy, 1) + size(obj.pxz, 1);
            
            % Initialize the resulting matrix R(X, Y, Z) to zeros
            obj.triplet_nonzero = zeros(obj.NX, obj.NY, obj.NZ);
            
            % % Create a 3D mask where we will set R(x, y, z) to 1 if both conditions are met
            % index = 1;
            % for i = 1:length(obj.nzXy)
            %     for j = 1:length(obj.nzXz)
            %         if obj.nzXy(i) == obj.nzXz(j)
            %             obj.triplet_nonzero(obj.nzXy(i), obj.nzY(i), obj.nzZ(j)) = index; %sub2ind(size(obj.triplet_nonzero), nzXy(i), nzY(i), nzZ(j));%1;
            %             index = index + 1;
            %         end
            %     end
            % end

            index = 1;
            for x=1:size(obj.b_xy,1)
                for y=1:size(obj.b_xy,2)
                    for z=1:size(obj.b_xz,2)
                        if (obj.b_xy(x,y)>0) && (obj.b_xz(x,z)>0)
                            obj.triplet_nonzero(x,y,z) = index;
                            index = index + 1;
                        end
                    end
                end
            end

            [obj.nzXYZx, obj.nzXYZyz, ~] = find(obj.triplet_nonzero);
            % obj.pxy = [nzXy, nzY, pxy];
            % obj.pxz = [nzXz, nzZ, pxz];

            % triplet_indices = ismember(obj.pxy(:,1), obj.pxz(:,1));
            % obj.triplet_nonzero = [obj.pxy(triplet_indices, :), obj.pxz(triplet_indices, 2:3)];
            % obj.n = size(obj.triplet_nonzero, 1);
            % obj.m = size(obj.pxy, 1) + size(obj.pxz, 1);
        end

        % function create_model(obj)
        %     % (c) Abdullah Makkeh, Dirk Oliver Theis
        %     % Permission to use and modify under Apache License version 2.0
        % 
        %     n_vars = 3 * obj.n;
        %     n_cons = obj.n + obj.m;
        % 
        %     obj.b = zeros(n_cons, 1,'double');
        % 
        %     Eqn = [];
        %     Var = [];
        %     Coeff = [];
        % 
        %     eqn = 0;
        %     for i = 1:obj.n
        %         x = obj.nzXYZx(i);
        %         [z, y] = ind2sub([obj.NZ, obj.NY],obj.nzXYZyz(i));
        %         % xyz = obj.triplet_nonzero(i,:);
        %         eqn = eqn + 1;
        %         p_var = obj.p_vidx(i);
        %         Eqn   = [Eqn;   eqn   ];
        %         Var   = [Var;   p_var ];
        %         Coeff = [Coeff; -1    ];
        % 
        %         % y = xyz(2);
        %         % z = xyz(4);
        % 
        %         subtriplet = obj.triplet_nonzero(:,y,z);
        %         select_idx = subtriplet(subtriplet>0);%find(obj.triplet_nonzero(:,y,z));%obj.triplet_nonzero(:, 2) == y & obj.triplet_nonzero(:, 4) == z;
        %         if ~isempty(select_idx)
        %             select_tab_qid = obj.q_vidx(select_idx);
        %             Eqn = [Eqn; repmat(eqn, length(select_tab_qid), 1)];
        %             Var = [Var; select_tab_qid(:)];
        %             Coeff = [Coeff; ones(length(select_tab_qid), 1)];
        %         end
        %     end
        % 
        %     for row = 1:length(obj.nzXy)%1:size(obj.pxy, 1)
        %         eqn = eqn + 1;
        %         x = obj.nzXy(row);
        %         y = obj.nzY(row);
        %         subtriplet = squeeze(obj.triplet_nonzero(x,y,:));
        %         selec_triplets = squeeze(subtriplet(subtriplet>0));
        %         % selec_triplets = triplet_nonzero(); %obj.triplet_nonzero(obj.triplet_nonzero(:, 1) == obj.pxy(row, 1) & obj.triplet_nonzero(:, 2) == obj.pxy(row, 2), :);
        %         if ~isempty(selec_triplets)
        %             selec_triplets_qid = obj.q_vidx(selec_triplets);
        %             Eqn = [Eqn; repmat(eqn, length(selec_triplets_qid), 1)];
        %             Var = [Var; selec_triplets_qid(:)];
        %             Coeff = [Coeff; ones(length(selec_triplets_qid), 1)];
        %         end
        % 
        %         obj.b(eqn) = obj.pxy(row);
        %     end
        % 
        %     for row = 1:length(obj.nzXz)
        %         eqn = eqn + 1;
        %         x = obj.nzXz(row);
        %         z = obj.nzZ(row);
        %         subtriplet = squeeze(obj.triplet_nonzero(x,:,z));
        %         selec_triplets = squeeze(subtriplet(subtriplet>0));
        %         if ~isempty(selec_triplets)
        %             selec_triplets_qid = obj.q_vidx(selec_triplets);
        %             Eqn = [Eqn; repmat(eqn, length(selec_triplets_qid), 1)];
        %             Var = [Var; selec_triplets_qid(:)];
        %             Coeff = [Coeff; ones(length(selec_triplets_qid), 1)];
        %         end
        % 
        %         obj.b(eqn) = obj.pxz(row);
        %     end
        % 
        %     obj.Eqn_list = Eqn;
        %     obj.A = sparse(Eqn, Var, Coeff, n_cons, n_vars);
        % 
        %     Ieq = [];
        %     IeqVar = [];
        %     Coeff = [];
        %     for i = 1:obj.n
        %         r_var = obj.r_vidx(i);
        %         q_var = obj.q_vidx(i);
        %         p_var = obj.p_vidx(i);
        % 
        %         Ieq = [Ieq; length(Ieq) + 1; length(Ieq) + 2; length(Ieq) + 3];
        %         IeqVar = [IeqVar; r_var; p_var; q_var];
        %         Coeff = [Coeff; -1; -1; -1];
        %     end
        % 
        %     obj.G = sparse(Ieq, IeqVar, Coeff, n_vars, n_vars);
        %     obj.h = zeros(n_vars, 1);
        %     obj.dims.e = obj.n;
        % 
        %     % Objective function:
        %     obj.c = zeros(n_vars, 1);
        %     for i = 1:obj.n
        %         obj.c(obj.r_vidx(i)) = -1;
        %     end
        % end

        function create_model(obj)
            % Initialize the number of variables and constraints
            n_vars = 3 * obj.n;  % 3 variables (r, p, q) per triplet
            n_cons = obj.n + obj.m;  % Number of constraints
        
            % Initialize the constraint vector 'b'
            obj.b = zeros(n_cons, 1, 'double');
        
            % Initialize the matrices to store the coefficients of the equations
            Eqn = [];  % Equation indices
            Var = [];  % Variable indices
            Coeff = [];  % Coefficient values
        
            eqn = 0;  % Equation counter
        
            % Loop through each non-zero triplet and create the equations
            for i = 1:obj.n
                % Get the indices of x, y, z for the current triplet
                % x = obj.nzXYZx(i);
                % [y, z] = ind2sub([obj.NY, obj.NZ], obj.nzXYZyz(i));
                [x, yz] = find(obj.triplet_nonzero==i);
                % [y, z] = ind2sub([obj.NY, obj.NZ], yz);
                y = mod(yz - 1, obj.NY) + 1;
                z = (yz - y) / obj.NY + 1;
                % Create the equation for p_i + q_j - r_k <= 0
                eqn = eqn + 1;
                p_var = obj.p_vidx(i);  % Get the variable index for p
                Eqn = [Eqn; eqn];
                Var = [Var; p_var];
                Coeff = [Coeff; -1];  % Coefficient for p
        
                % Find the corresponding q variables
                subtriplet = obj.triplet_nonzero(:, y, z);  % Select the subtriplet with x fixed
                select_idx = subtriplet(subtriplet > 0);  % Find the non-zero entries
        
                if ~isempty(select_idx)
                    select_tab_qid = obj.q_vidx(select_idx);  % Get the q indices
                    Eqn = [Eqn; repmat(eqn, length(select_tab_qid), 1)];
                    Var = [Var; select_tab_qid(:)];
                    Coeff = [Coeff; ones(length(select_tab_qid), 1)];  % Coefficients for q
                end
            end
        
            % Add the constraints for the marginal distributions p(x,y)
            for row = 1:length(obj.nzXy)
                eqn = eqn + 1;
                x = obj.nzXy(row);
                y = obj.nzY(row);
        
                % Get the subtriplet for the current (x, y) pair
                % subtriplet = squeeze(obj.triplet_nonzero(x, y, :));
                % selec_triplets = squeeze(subtriplet(subtriplet > 0));
                subtriplet = obj.triplet_nonzero(x, y, :);
                selec_triplets = subtriplet(subtriplet > 0);
        
                if ~isempty(selec_triplets)
                    selec_triplets_qid = obj.q_vidx(selec_triplets);
                    Eqn = [Eqn; repmat(eqn, length(selec_triplets_qid), 1)];
                    Var = [Var; selec_triplets_qid(:)];
                    Coeff = [Coeff; ones(length(selec_triplets_qid), 1)];
                end
        
                obj.b(eqn) = obj.pxy(row);  % Assign the probability to the corresponding equation
            end
        
            % Add the constraints for the marginal distributions p(x,z)
            for row = 1:length(obj.nzXz)
                eqn = eqn + 1;
                x = obj.nzXz(row);
                z = obj.nzZ(row);
        
                % Get the subtriplet for the current (x, z) pair
                subtriplet = squeeze(obj.triplet_nonzero(x, :, z));
                selec_triplets = squeeze(subtriplet(subtriplet > 0));
        
                if ~isempty(selec_triplets)
                    selec_triplets_qid = obj.q_vidx(selec_triplets);
                    Eqn = [Eqn; repmat(eqn, length(selec_triplets_qid), 1)];
                    Var = [Var; selec_triplets_qid(:)];
                    Coeff = [Coeff; ones(length(selec_triplets_qid), 1)];
                end
        
                obj.b(eqn) = obj.pxz(row);  % Assign the probability to the corresponding equation
            end
        
            % Convert the collected data into a sparse matrix 'A'
            obj.Eqn_list = Eqn;
            obj.Var_list = Var;
            obj.Coeff_list = Coeff;
            obj.A = sparse(Eqn, Var, Coeff, n_cons, n_vars);
        
            % Set up the inequality constraints (G * rpq <= h)
            Ieq = [];
            IeqVar = [];
            Coeff = [];
            for i = 1:obj.n
                r_var = obj.r_vidx(i);
                q_var = obj.q_vidx(i);
                p_var = obj.p_vidx(i);
        
                Ieq = [Ieq; length(Ieq) + 1; length(Ieq) + 2; length(Ieq) + 3];
                IeqVar = [IeqVar; r_var; p_var; q_var];
                Coeff = [Coeff; -1; -1; -1];  % Set up the coefficients for the inequalities
            end
        
            % Convert to sparse matrix 'G'
            obj.G = sparse(Ieq, IeqVar, Coeff, n_vars, n_vars);
            obj.h = zeros(n_vars, 1);  % Initialize 'h' as a zero vector
            obj.dims.e = obj.n;  % Set the dimensions for the exponential cone constraints
        
            % Objective function (minimizing the sum of r_i)
            obj.c = zeros(n_vars, 1);
            for i = 1:obj.n
                obj.c(obj.r_vidx(i)) = -1;  % Objective function minimizes 'r'
            end
        end

        
        function v_idx = r_vidx(obj, i)
            v_idx = 3 * (i - 1) + 1;
        end

        function v_idx = p_vidx(obj, i)
            v_idx = 3 * (i - 1) + 2;
        end

        function v_idx = q_vidx(obj, i)
            v_idx = 3 * (i - 1) + 3;
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
            [obj.sol_rpq, obj.sol_lambda, obj.sol_info, obj.sol_slack, obj.sol_mu] = ecos(obj.c, obj.G, obj.h, obj.dims, obj.A, obj.b, obj.ecos_kwargs);
            
            if size(obj.sol_rpq,1)> 0
                obj.solution.r = obj.sol_rpq(1:3:end);
                obj.solution.p = obj.sol_rpq(2:3:end);
                obj.solution.q = obj.sol_rpq(3:3:end);
                result = 'success';
            else
                % "x" not in struct solution
                result = 'x not in struct solution -- No Solution Found!!!';
            end
        end

        function result = entropy_X2(obj)
            px = sum(obj.b_xy,2);
            result = - dot(nonzeros(px), log(nonzeros(px)));
        end

        function addqdistr(obj)
            obj.triplet_nonzeroqid = 3 * obj.triplet_nonzero(obj.triplet_nonzero>0);
            obj.triplet_nonzeroq = zeros(size(obj.triplet_nonzero)); %obj.sol_rpq(obj.triplet_nonzeroqid);
            % index=1;
            % for x=1:size(obj.b_xy,1)
            %     for y=1:size(obj.b_xy,2)
            %         for z=1:size(obj.b_xz,2)
            %             if (obj.b_xy(x,y)>0) && (obj.b_xz(x,z)>0)
            %                 obj.triplet_nonzeroq(x,y,z) = obj.solution.q(index);
            %                 index = index + 1;
            %             end
            %         end
            %     end
            % end
            for i = 1:obj.n
                [x, yz] = find(obj.triplet_nonzero==i);
                [y, z] = ind2sub([obj.NY, obj.NZ], yz);
                obj.triplet_nonzeroq(x, y, z) = obj.solution.q(i);
            end
        end

        function create_qxyz(obj)
            obj.triplet_nonzeroqpos = obj.triplet_nonzeroq;
            obj.triplet_nonzeroqpos(obj.triplet_nonzeroqpos<0) = 0;
            
            obj.qyz = sum(obj.triplet_nonzeroqpos, 1);
            obj.qyz = repmat(obj.qyz,obj.NX, 1,1);
            %qy = squeeze(sum(obj.qyz,2));
            %qz = squeeze(sum(obj.qyz,1));
        end

        function result=condentropy2(obj)
            mask = obj.triplet_nonzeroqpos > 0;
            condentr = -obj.triplet_nonzeroqpos(mask) .* log(obj.triplet_nonzeroqpos(mask) ./ obj.qyz(mask));
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
            result = sum(tmp_table{:,'condentrorig'},1);
        end
        
        function result = condentropy_orig2(obj, pdf)
            pyz = sum(pdf, 1);
            pyz = repmat(pyz,obj.NX, 1, 1);
            tmp_condentrorig = -pdf(pdf>0) .* log(pdf(pdf>0) ./ pyz(pdf>0));
            result = sum(tmp_condentrorig,'all');
        end
        

        function result = condYmutinf2(obj)
            mask = obj.triplet_nonzeroqpos > 0;
            obj.qy = sum(obj.triplet_nonzeroqpos, [1 3]);
            obj.qy = repmat(obj.qy, obj.NX, 1, obj.NZ);
            b_xy_repeated = repmat(obj.b_xy, [1 1 obj.NZ]);
            condYmutinf = obj.triplet_nonzeroqpos(mask) .* log(obj.triplet_nonzeroqpos(mask) .* obj.qy(mask) ./ (b_xy_repeated(mask) .* obj.qyz(mask)));
            result = sum(condYmutinf);
        end

        function result = condZmutinf2(obj)
            mask = obj.triplet_nonzeroqpos > 0;
            obj.qz = sum(obj.triplet_nonzeroqpos, [1 2]);
            obj.qz = repmat(obj.qz, obj.NX, obj.NY, 1);
            b_xz_repeated = repmat(obj.b_xz, obj.NY);
            condZmutinf = obj.triplet_nonzeroqpos(mask) .* log(obj.triplet_nonzeroqpos(mask) .* obj.qz(mask) ./ (b_xz_repeated(mask) .* obj.qyz(mask)));
            result = sum(condZmutinf);
        end

        function result = condZmutinf3(obj)
            MI = 0;

            % Loop through the non-zero entries in the triplet matrix
            for i = 1:obj.n
                [x, yz] = find(obj.triplet_nonzero==i);
                [y, z] = ind2sub([obj.NY, obj.NZ], yz);
        
                % Calculate the joint probability p(x, y, z)
                p_xyz = obj.triplet_nonzeroqpos(x, y, z);
        
                if p_xyz > 0
                    % Calculate the marginal probabilities
                    p_yz = obj.qyz(1, y, z); % Marginal probability p(y, z)
                    p_xz = obj.pxz(obj.nzXz == x & obj.nzZ == z); % Marginal probability p(x, z)
                    p_z = sum(obj.pxz(obj.nzZ == z)); % Marginal probability p(z)
        
                    % Mutual information contribution for the current triplet
                    MI = MI + p_xyz * log(p_xyz * p_z / (p_xz * p_yz));
                end
            end
        
            % Result is the mutual information
            result = MI;
        end
    end
    
end
