classdef pid_lattice < handle
    % PID_LATTICE Constructs and operates on a redundancy lattice for probabilistic variables.
    %
    % This class builds and evaluates a redundancy lattice using Imin or IMMI/MMI
    % as the redundancy measure. If a precomputed lattice file is not found,
    % the lattice is generated dynamically.
    %
    % It supports:
    %   • Discrete (non-Gaussian) case: work with a full joint probability mass
    %     function (pdf) over [sources, target].
    %   • Gaussian case: work with a covariance matrix and compute mutual
    %     information using Gaussian closed forms.
    %
    % -------------------------------------------------------------------------
    % Properties
    % -------------------------------------------------------------------------
    %   lat                - 2×N cell array storing lattice nodes and (optionally) values.
    %                        lat(1, :) are nodes; each node is a cell array of source-index arrays.
    %   pi                 - 1×N vector of PID atom values (partial informations) per node.
    %   pdf                - Discrete: joint probability tensor over sources and target.
    %                        Gaussian: covariance matrix over [all sources, target].
    %   nsources           - Number of source variables.
    %   is_gaussian        - If true, use Gaussian mutual information (for IMMI/MMI path).
    %   redundancy_measure - 'Imin' or 'IMMI' (default: 'Imin').
    %   source_dims        - (Gaussian) 1×nsources cell, each entry contains indices for a source.
    %   target_dims        - (Gaussian) Index/indices for the target variable.
    %
    % -------------------------------------------------------------------------
    % Constructor
    % -------------------------------------------------------------------------
    %   pid_lattice(nsources, redundancy_measure, is_gaussian)
    %     Initializes the object and attempts to load a cached lattice file
    %     ('lat%dSources.mat'); if not found, generates all valid redundancy
    %     lattice nodes.
    %
    % -------------------------------------------------------------------------
    % Lattice utilities
    % -------------------------------------------------------------------------
    %   power_set(var_list)
    %     Return the (non-empty) power set of var_list as a cell array of index arrays.
    %
    %   is_node_red(x)
    %     Return true if the candidate node x (a cell of source-index arrays)
    %     is valid for the redundancy lattice (i.e., no element is a subset of another).
    %
    %   node_issubsetany(x, i, j)
    %     Return true if x{i} is not a subset/superset of x{j} (or i==j).
    %
    %   node_issubset(x, y)
    %     Return true if every element in node x is a subset of at least one
    %     element in node y (partial order on redundancy lattice).
    %
    %   node_issame(x, y)
    %     Return true if nodes x and y contain the same elements (order-insensitive).
    %
    %   get_down(node)
    %     Return all nodes in the lattice that are ≤ node in the partial order
    %     (including the node itself).
    %
    %   get_strict_down(node)
    %     Return all nodes strictly below the given node in the partial order.
    %
    % -------------------------------------------------------------------------
    % PID / redundancy evaluation
    % -------------------------------------------------------------------------
    %   calculate_latvals(p_distr)
    %     Compute PID atom values (obj.pi) for all lattice nodes using either
    %     Imin (discrete) or IMMI/MMI (Gaussian or if redundancy_measure == 'IMMI').
    %     Returns obj.pi.
    %
    %   calculate_atom(p_distr, target, sources)
    %     Recursive Möbius inversion on the lattice for the discrete case.
    %     Base redundancy term uses Imin (if redundancy_measure == 'Imin') or
    %     MMI (if redundancy_measure == 'IMMI'); subtract strict-down atoms.
    %
    % -------------------------------------------------------------------------
    % Redundancy measures
    % -------------------------------------------------------------------------
    %   Imin(p_distr, target, sources)
    %     Minimum specific information between target and the set of sources,
    %     computed from the discrete joint pdf. Requires helper function
    %     create_prob_ts to build pairwise joints.
    %
    %   MMI(p_distr, target, sources)
    %     Minimum of mutual informations I(S_i; T) over all source collections
    %     in 'sources'. Uses:
    %       • Discrete: marginalize() + mutualInformationLast().
    %       • Gaussian: gaussianMI() with source/target index sets.
    %
    %   specific_information(p_distr, specific_val_dim, specific_val_index)
    %     Specific information I(X=x; Y) for a fixed value index along dimension
    %     specific_val_dim of the joint pdf p_distr (log base 2).
    %
    % -------------------------------------------------------------------------
    % Information-theoretic helpers
    % -------------------------------------------------------------------------
    %   mutualInformationLast(p)
    %     Compute I(X; Z) where the last dimension of tensor p is Z and all
    %     preceding dimensions comprise X. Returns MI in bits.
    %
    %   marginalize(p, keep_dims)
    %     Sum out all dimensions of p not listed in keep_dims. Returns the
    %     squeezed marginal tensor.
    %
    %   gaussianMI(cov_mat, source_idx, target_idx)
    %     Compute Gaussian mutual information I(X; Y) = 0.5*log2(|Σ_X||Σ_Y|/|Σ_{XY}|),
    %     where indices select blocks from cov_mat. Small diagonal jitter is added
    %     for numerical stability.
    %
    % Notes:
    %   • Discrete vs Gaussian:
    %       - Discrete workflows assume p_distr is a normalized joint pdf tensor.
    %       - Gaussian workflows assume p_distr is a covariance matrix with
    %         variable ordering [sources, target].
    %   • Lattice nodes are expressed as cells of source-index arrays; e.g.,
    %     {{1}, {2,3}} represents the node with two “inputs”: source 1 and the
    %     joint source {2,3}.
    %

    % Copyright (C) 2024 Gabriel Matias Lorenz, Nicola Marie Engel
    % This file is part of MINT.
    % This program is free software: you can redistribute it and/or modify
    % it under the terms of the GNU General Public License as published by
    % the Free Software Foundation, either version 3 of the License, or
    % (at your option) any later version.
    %
    % This program is distributed in the hope that it will be useful,
    % but WITHOUT ANY WARRANTY; without even the implied warranty of
    % MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    % GNU General Public License for more details.
    %
    % You should have received a copy of the GNU General Public License
    % along with this program.  If not, see <http://www.gnu.org/licenses

    properties
        lat
        pi
        pdf
        nsources
        is_gaussian = false
        redundancy_measure = 'Imin' % 'Imin' or 'IMMI'
        source_dims     % Cell array of source variable indices (for Gaussian)
        target_dims
    end

    methods
        function obj = pid_lattice(nsources, redundancy_measure, is_gaussian)
            obj.nsources = nsources;
            if nargin >= 2
                obj.redundancy_measure = redundancy_measure;
            end
            if nargin >= 3
                obj.is_gaussian = is_gaussian;
            end

            if obj.is_gaussian
                obj.source_dims = arrayfun(@(s) {s}, 1:obj.nsources);  % Each source is 1D by default
                obj.target_dims = obj.nsources + 1;                    % Target index
            end


            filename = sprintf('lat%dsources', nsources);
            try
                lat = load(join([filename '.mat']));
                obj.lat = lat.(filename);
            catch
                fprintf('[INFO] Lattice file not found. Generating lattice...\n');
                p_set = obj.power_set(1:nsources);
                pp_set = obj.power_set(p_set);
                valid_nodes = {};
                for k = 1:numel(pp_set)
                    node = pp_set{k};
                    if obj.is_node_red(node)
                        valid_nodes{end+1} = node;
                    end
                end
                obj.lat = cell(2, numel(valid_nodes));
                for i = 1:numel(valid_nodes)
                    obj.lat{1, i} = valid_nodes{i};
                end
                fprintf('[INFO] Lattice generated with %d nodes.\n', numel(valid_nodes));
            end
        end

        function result = power_set(~, var_list)
            n = length(var_list);
            result = cell(1, 2^n - 1);
            for i = 1:(2^n - 1)
                result{i} = var_list(bitget(i, 1:n) == 1);
            end
        end

        function validity = is_node_red(obj, x)
            validity = true;
            flat = [x{:}];
            if any(flat > obj.nsources)
                validity = false; return;
            end
            for i = 1:numel(x)
                for j = 1:numel(x)
                    if i ~= j && (all(ismember(x{i}, x{j})) || all(ismember(x{j}, x{i})))
                        validity = false;
                        return;
                    end
                end
            end
        end

        function issubsetany = node_issubsetany(obj, x, i, j)
            if i~=j
                issubsetany = ~(all(ismember(x{i}, x{j})) | all(ismember(x{j}, x{i})));
            else
                issubsetany = true;
            end
        end

        function inclusion = node_issubset(obj, x, y)
            % The function is_min_red determines if a node x precedes another node y
            % according to the partial order in a redundancy lattice. It checks whether
            % every element in node x is a subset of at least one element in node y.
            %
            % Input:
            %     x: A cell array representing a node in the redundancy lattice.
            %     y: A cell array representing another node in the redundancy lattice.

            % Output:
            %    inclusion: A logical value indicating whether every element in node
            %               x is a subset of at least one element in node y
            %               (true if included, false otherwise).
            inclusion = true;
            yi=1;
            while inclusion==true && yi<=numel(y)
                node_inclusion = arrayfun(@(xi) all(ismember(x{xi}, y{yi})),1:numel(x));
                % node_inclusion = arrayfun(@(xi) all(ismember(x{xi}{1}, y{yi}{1})),1:numel(x));
                yi=yi+1;
                inclusion = any(node_inclusion);
            end
        end

        function down_list = get_down(obj, node)
            down_mask = arrayfun(@(j) obj.node_issubset(obj.lat{1,j}, node), 1:size(obj.lat,2), 'UniformOutput', true);
            % down_mask = arrayfun(@(j) obj.node_issubset(obj.lat(1,j), node), 1:size(obj.lat,2), 'UniformOutput', true);
            down_list = obj.lat(1,down_mask);
        end

        function strict_down_list = get_strict_down(obj, node)
            down_list = obj.get_down(node);
            strict_down_mask = arrayfun(@(j) ~obj.node_issame(down_list{j}, node), 1:numel(down_list), 'UniformOutput', true);
            % strict_down_mask = arrayfun(@(j) ~obj.node_issame(down_list{j}, node{1}), 1:numel(down_list), 'UniformOutput', true);
            strict_down_list = down_list(strict_down_mask);
        end

        function equality = node_issame(~, x, y)
            if numel(x) ~= numel(y)
                equality = false;
            else
                equality = all(arrayfun(@(i) isequal(sort(x{i}), sort(y{i})), 1:numel(x)));
            end
        end

        function latvals = calculate_latvals(obj, p_distr)
            obj.pdf = p_distr;
            obj.pi = zeros(1, size(obj.lat,2));
            for node = 1:size(obj.lat,2)
                sources = obj.lat{1,node};
                target = obj.nsources + 1;
                if obj.is_gaussian
                    obj.pi(node) = obj.calculate_atom_mmi(p_distr, target, sources);
                else
                    obj.pi(node) = obj.calculate_atom(p_distr, target, sources);
                end
            end
            latvals = obj.pi;
        end

        function atom = calculate_atom(obj, p_distr, target, sources)
            if obj.is_gaussian | strcmp(obj.redundancy_measure, 'IMMI')
                atom = obj.Imin(p_distr, target, sources);
            else
                atom = obj.MMI(p_distr, target, sources);
            end
            down_nodes = obj.get_strict_down(sources);
            % atom = imin;
            for i = 1:numel(down_nodes)
                atom = atom - obj.calculate_atom(p_distr, target, down_nodes{i});
            end
        end

        function val = Imin(obj, p_distr, target, sources)
            if exist('create_prob_ts', 'file') ~= 2
                error('Missing function: create_prob_ts must be defined or on path.');
            end
            spec_info = zeros(size(p_distr, target), numel(sources));
            for s_i = 1:numel(sources)
                p_joint = create_prob_ts(p_distr, [target, sources{s_i}]);
                for t_i = 1:size(p_joint,1)
                    spec_info(t_i, s_i) = obj.specific_information(p_joint, 1, t_i);
                end
            end
            p_target = squeeze(sum(p_distr, setdiff(1:ndims(p_distr), target)));
            val = p_target' * min(spec_info, [], 2);
        end

        function mmi_val = MMI(obj, p_distr, target, sources)
            mmi_vals = zeros(1, numel(sources));
            for i = 1:numel(sources)
                if obj.is_gaussian
                    source_idx = cell2mat(obj.source_dims(sources{i}));
                    target_idx = obj.target_dims;
                    mmi_vals(i) = obj.gaussianMI(p_distr, source_idx, target_idx);
                else
                    marg = obj.marginalize(p_distr, [sources{i}, target]);
                    mmi_vals(i) = obj.mutualInformationLast(marg);
                end
            end
            mmi_val = min(mmi_vals);
        end
        
        function specinf = specific_information(~, p_distr, dim, val_index)
            bits = 1/log(2);
            permuted = permute(p_distr, [dim setdiff(1:ndims(p_distr), dim)]);
            marg_dim = sum(permuted, 2:ndims(p_distr));
            p_val = marg_dim(val_index);
            marg_other = squeeze(sum(permuted, 1));
            permuted = reshape(permuted, size(p_distr, dim), []);
            specinf = 0;
            for i = 1:size(permuted, 2)
                if permuted(val_index, i) > 0
                    quot = permuted(val_index, i) / p_val;
                    specinf = specinf + quot * log(quot / marg_other(i));
                end
            end
            specinf = specinf * bits;
        end

        function MI = mutualInformationLast(~, p)
            p = p / sum(p(:));
            z_dim = ndims(p);
            x_dims = 1:(z_dim - 1);
            pz = sum(p, x_dims);
            px = sum(p, z_dim);
            px_exp = repmat(px, [ones(1, z_dim - 1), size(p, z_dim)]);
            pz_exp = repmat(reshape(pz, [ones(1, z_dim - 1), numel(pz)]), size(px));
            valid = p > 0 & px_exp > 0 & pz_exp > 0;
            MI = sum(p(valid) .* log2(p(valid) ./ (px_exp(valid) .* pz_exp(valid))));
        end

        function marginal = marginalize(~, p, keep_dims)
            all_dims = 1:ndims(p);
            sum_dims = setdiff(all_dims, keep_dims);
            for d = sum_dims
                p = sum(p, d);
            end
            marginal = squeeze(p); %permute(p, sort(keep_dims));
        end


        function MI = gaussianMI(~, cov_mat, source_idx, target_idx)
            all_idx = [source_idx, target_idx];
            cov_xy = cov_mat(all_idx, all_idx);
            cov_x = cov_mat(source_idx, source_idx);
            cov_y = cov_mat(target_idx, target_idx);
        
            cov_x = cov_x + 1e-10 * eye(size(cov_x));
            cov_y = cov_y + 1e-10 * eye(size(cov_y));
            cov_xy = cov_xy + 1e-10 * eye(size(cov_xy));
        
            MI = 0.5 * log2(det(cov_x) * det(cov_y) / det(cov_xy));
        end


    end
end
