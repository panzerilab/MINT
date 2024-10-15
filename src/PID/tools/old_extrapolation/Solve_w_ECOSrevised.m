classdef Solve_w_ECOSrevised < handle
    properties
        ecos_kwargs         
        verbose             
        c                   
        G
        h
        dims                
        A                   
        b
        sol_rpq             
        sol_slack
        sol_lambda
        sol_mu
        sol_info
        solution            
        b_xy                
        b_xz
        NX                  
        NY
        NZ
        triplet_nonzero     
        n                   
        m
        table_pxy           
        table_pxz
        Eqn_list            
    end

    methods
        function obj = Solve_w_ECOSrevised(marg_xy, marg_xz)
            obj.ecos_kwargs = struct();
            obj.verbose = false;

            obj.c = [];
            obj.G = [];
            obj.h = [];
            obj.dims = struct();
            obj.A = [];
            obj.b = [];

            obj.sol_rpq = [];
            obj.sol_slack = [];
            obj.sol_lambda = [];
            obj.sol_mu = [];
            obj.sol_info = [];
            obj.solution= struct;

            obj.b_xy = marg_xy;
            obj.b_xz = marg_xz;
            obj.NX = size(obj.b_xy, 1);
            obj.NY = size(obj.b_xy, 2);
            obj.NZ = size(obj.b_xz, 2);

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

            obj.table_pxy = array2table([nzXy nzY pxy], 'VariableNames',{'X', 'Y', 'pxy'} );
            obj.table_pxz = array2table([nzXz nzZ pxz], 'VariableNames',{'X', 'Z', 'pxz'} );

            [~, idx_pxy, idx_pxz] = intersect([nzXy nzY], [nzXz nzZ], 'rows');
            obj.triplet_nonzero = [nzXy(idx_pxy) nzY(idx_pxy) nzZ(idx_pxz)];
            % Convert the key columns from both tables to arrays
            keys1 = [obj.table_pxy.X];
            keys2 = [obj.table_pxz.X];
            
            % Find the matching rows using ismember
            [isCommon, idxTable1, idxTable2] = ismember(keys1, keys2, 'rows');
            
            % Select the rows that are common in both tables
            result_table1 = obj.table_pxy(idxTable1(isCommon), :);
            result_table2 = obj.table_pxz(idxTable2(isCommon), :);
            
            % Combine the results if needed (e.g., concatenate columns)
            result = [result_table1, result_table2];
            obj.triplet_nonzero(:, 4) = (1:size(obj.triplet_nonzero, 1))';
            obj.n = size(obj.triplet_nonzero, 1);
            obj.m = size(obj.table_pxy, 1) + size(obj.table_pxz, 1);
        end

        function create_model(obj)
            n_vars = 3 * obj.n;
            n_cons = obj.n + obj.m;

            obj.b = zeros(n_cons, 1,'double');
            Eqn = [];
            Var = [];
            Coeff = [];

            eqn = 0;
            for i = 1:obj.n
                eqn = eqn+1;
                p_var = obj.p_vidx(i);
                Eqn = [Eqn; eqn];
                Var = [Var; p_var];
                Coeff = [Coeff; -1];

                y = obj.triplet_nonzero(i, 2);
                z = obj.triplet_nonzero(i, 3);

                select_idx = obj.triplet_nonzero(:, 2) == y & obj.triplet_nonzero(:, 3) == z;
                if any(select_idx)
                    q_vars = obj.q_vidx(obj.triplet_nonzero(select_idx, 4));
                    Eqn = [Eqn; repmat(eqn, sum(select_idx), 1)];
                    Var = [Var; q_vars];
                    Coeff = [Coeff; ones(sum(select_idx), 1)];
                end
            end

            for row = 1:size(obj.table_pxy, 1)
                eqn = eqn + 1;
                select_idx = ismember(obj.triplet_nonzero(:, 1:2), obj.table_pxy(row, 1:2), 'rows');
                if any(select_idx)
                    q_vars = obj.q_vidx(obj.triplet_nonzero(select_idx, 4));
                    Eqn = [Eqn; repmat(eqn, sum(select_idx), 1)];
                    Var = [Var; q_vars];
                    Coeff = [Coeff; ones(sum(select_idx), 1)];
                end
                obj.b(eqn) = obj.b_xy(obj.table_pxy(row, 1), obj.table_pxy(row, 2));
            end

            for row = 1:size(obj.table_pxz, 1)
                eqn = eqn + 1;
                select_idx = ismember(obj.triplet_nonzero(:, [1, 3]), obj.table_pxz(row, [1, 2]), 'rows');
                if any(select_idx)
                    q_vars = obj.q_vidx(obj.triplet_nonzero(select_idx, 4));
                    Eqn = [Eqn; repmat(eqn, sum(select_idx), 1)];
                    Var = [Var; q_vars];
                    Coeff = [Coeff; ones(sum(select_idx), 1)];
                end
                obj.b(eqn) = obj.b_xz(obj.table_pxz(row, 1), obj.table_pxz(row, 2));
            end

            obj.Eqn_list = Eqn;
            obj.A = sparse(Eqn, Var, Coeff, n_cons, n_vars);

            Ieq = [];
            IeqVar = [];
            Coeff = [];
            for i = 1:obj.n
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
            if ~isempty(obj.verbose)
                obj.ecos_kwargs.verbose = obj.verbose;
            end
            obj.ecos_kwargs.feastol = 1E-7;
            obj.ecos_kwargs.abstol = 1E-6;
            obj.ecos_kwargs.reltol = 1E-6;

            [x, y, info, s, z] = ecos(obj.c, obj.G, obj.h, obj.dims, obj.A, obj.b, obj.ecos_kwargs);

            if size(x, 1) > 0
                obj.sol_rpq = x;
                obj.sol_slack = s;
                obj.sol_lambda = y;
                obj.sol_mu = z;
                obj.sol_info = info;
                result = 'success';
            else
                result = 'No Solution Found';
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
            [gx, ~] = findgroups(obj.table_pxy(:, 1));
            p_x = splitapply(@sum, obj.table_pxy(:, 3), gx);
            result = -sum(p_x .* log(p_x));
            % obj.triplet_nonzero(:, 5) = p_x;
        end

        function result = condentropy(obj)
            result = 0;
            for yv = 1:obj.NY
                for zv = 1:obj.NZ
                    select_idx = obj.triplet_nonzero(:, 2) == yv & obj.triplet_nonzero(:, 3) == zv;
                    marg_x = sum(max(0, obj.sol_rpq(obj.q_vidx(obj.triplet_nonzero(select_idx, 4)))));
                    if marg_x > 0
                        q_vals = obj.sol_rpq(obj.q_vidx(obj.triplet_nonzero(select_idx, 4)));
                        result = result - sum(q_vals .* log(q_vals / marg_x));
                    end
                end
            end
        end

        function addqdistr(obj)
            obj.triplet_nonzero(:, 6) = obj.sol_rpq(obj.q_vidx(obj.triplet_nonzero(:, 4)));
        end

        function create_qxyz(obj)
            obj.triplet_nonzero(:, 7) = max(0, obj.triplet_nonzero(:, 6));

            gyz = accumarray(obj.triplet_nonzero(:, [2, 3]), obj.triplet_nonzero(:, 7));
            gy = accumarray(obj.triplet_nonzero(:, 2), obj.triplet_nonzero(:, 7));
            gz = accumarray(obj.triplet_nonzero(:, 3), obj.triplet_nonzero(:, 7));

            % qxyz = reshape(obj.triplet_nonzero(:, 7), [obj.NX, obj.NY, obj.NZ]);
            obj.solution.qxyz = obj.triplet_nonzero(:, 7);
            obj.solution.qyz = gyz;
            obj.solution.qy = gy;
            obj.solution.qz = gz;
        end
    end
end
