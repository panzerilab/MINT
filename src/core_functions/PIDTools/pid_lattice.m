classdef pid_lattice < handle 
% Description:
% This class implements a redundancy lattice, a mathematical structure used to model redundancy relationships between sets of variables. It is particularly designed for use in the context of probabilistic information theory.
% 
% Properties:
% - lat: A cell array representing the lattice nodes.
% - red: Redundancy values associated with lattice nodes.
% - pi: Probability distribution array.
% - pdf: Probability density function associated with the lattice nodes.
% 
% Methods:
% Constructor: pid_lattice(nsources)
%    - Initializes the pid_lattice object with an empty lattice and generates lattice nodes based on the given number of variables (nsources).
% 
% power_set(obj, var_list)
%    - Computes the power set of a given set of variables using bitwise operations in MATLAB.
% 
% is_node_red(obj, x)
%    - Checks if a given collection of nodes x represents a valid node in the redundancy lattice.
% 
% node_issubsetany(obj, x, i, j)
%    - Checks if node x[i] is a subset of any other node x[j] within a given collection.
% 
% node_issubset(obj, x, y)
%    - Determines if a node x precedes another node y in the redundancy lattice.
% 
% node_issame(obj, x, y)
%    - Checks if two nodes x and y are identical in terms of their elements.
% 
% get_down(obj, node)
%    - Returns the set of lattice nodes that are below the given node in the lattice.
% 
% get_strict_down(obj, node)
%    - Returns the set of lattice nodes that are strictly below the given node in the lattice.
% 
% get_red(obj, node)
%    - Computes the redundancy value associated with the given node.
% 
% Imin(obj, p_distr, target, sources)
%     - Computes the minimum specific information between a target variable and multiple source variables, given a probability distribution.
% 
% specific_information(obj, p_distr, specific_val_dim, specific_val_index)
%     - Calculates specific information between two variables based on a probability distribution.
% 
% Note: The class employs various methods to handle lattice nodes, check redundancy relationships, and calculate information-theoretic measures.

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
        red
        pi
        pdf
        nsources
    end

    methods
        function obj = pid_lattice(nsources)
            obj.nsources = nsources;
            filename= sprintf('lat%dsources',nsources);
            lat = load(join([filename '.mat']));
            obj.lat = lat.(filename);
            
            % obj.lat = {};
            % p_set = obj.power_set(1:nsources);
            % % p_set = p_set(2:end);
            % pp_set = obj.power_set(p_set);
            % % pp_set = pp_set(2:end);
            % lat_nodes = arrayfun(@(k) is_node_red(pp_set{k}), 1:numel(pp_set), 'UniformOutput', true);
            % obj.lat = pp_set(lat_nodes);
            % for i=1:size(obj.lat,2)
            %     obj.lat{2,i} = i;
            % end
            % 
        end

        function result = power_set(obj, var_list)
            % Power set using bitwise operations in MATLAB
            n = length(var_list);
            num_subsets = 2^n;

            % Preallocate result cell array
            result = cell(1, num_subsets - 1);

            % Loop through all possible binary representations
            for i = 1:(num_subsets - 1)
                result{i} = var_list(bitget(i, 1:n)==1);  % Extract subset
            end
        end

        function validity = is_node_red(obj, x)
            % The function is_node_red checks if a given collection of nodes x represents
            % a valid node in the redundancy lattice. A node is considered valid if none of
            % its elements are subsets of any other elements within the collection.
            %
            % Input:
            %     x: A cell array representing a collection of nodes in the redundancy lattice.
            %
            % Output:
            %     validity: A logical value indicating whether the input collection x represents a
            %               valid node in the redundancy lattice (true for valid, false otherwise).
           
            validity = true;
           
            x_tmp = [x{:}];
            if any(x_tmp > obj.nsources)
                validity = false;
                return
            end 
            
            i = 1;
            while validity && i <= numel(x)
                source_exclusion = arrayfun(@(j) obj.node_issubsetany(x, i, j), 1:numel(x), 'UniformOutput', true);
                validity = all(source_exclusion);
                i = i + 1;
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

        function equality = node_issame(obj, x, y)
            if numel(x) ==numel(y)
                equality = all(arrayfun(@(xi) all(ismember(x{xi}, y{xi})) && all(ismember(y{xi}, x{xi})),1:numel(x)));
            else
                equality = false;
            end
        end

        function get_red(obj)
            if isempty(obj.pdf)
                warning('You need to define a probability distribution function for your object.')
                return
            end 

            for i=1:size(obj.lat,2)
                obj.lat{2,i} = obj.Imin(obj.pdf, 1, obj.lat(1,i));
            end
        end

        function get_pi(obj, pdf)
            obj.pdf = pdf;
            obj.get_red();
            for i=1:size(obj.lat,2)
                strict_down = obj.get_strict_down(obj.lat(1,i));
                obj.lat{3,i} = obj.lat{2,i};
                for j=1:size(obj.lat,2)
                    if i==j
                        continue
                    end
                    if ismember(obj.lat{1,j}, strict_down)
                        obj.lat{3,i} = obj.lat{3,i} - obj.lat{2,j};
                    end
                end
            end
        end


        function imin_v = Imin(obj, p_distr, target, sources)
            % The function Imin computes the minimum specific information between a target variable
            % and multiple source variables given a probability distribution. Specific information
            % measures the information content of one variable with respect to another. The function
            % uses the specific_information and create_prob_ts functions to calculate specific information
            % for different source-target pairs.
            %
            % Inputs:
            %     p_distr: A probability distribution represented as a multi-dimensional array.
            %     target: The index representing the target variable within the probability distribution.
            %     sources: A cell array containing indices or names of source variables.
            %
            % Output:
            %     imin_v: The minimum specific information between the target variable and the sources.

            spec_inf_array = zeros(size(p_distr,target),numel(sources));
            for s_i = 1:numel(sources)
                % p_ts = create_prob_ts(p_distr, target, sources{s_i}{1});
                % p_ts = create_prob_ts(p_distr, target, sources{s_i});
                p_ts = create_prob_ts(p_distr, sources{s_i});
                for t_i = 1:size(p_ts,1)
                    spec_inf_array(t_i, s_i) = obj.specific_information(p_ts, 1, t_i);
                end
            end
            % vt = 1:ndims(p_distr);
            % vt([1, target]) = vt([target, 1]);
            % permuted_prob = permute(p_distr, vt);
            p_target = squeeze(sum(p_distr,setdiff(1:ndims(p_distr), target))); %sum(permuted_prob, 2:ndims(p_distr));
            imin_v = p_target' * squeeze(min(spec_inf_array, [], 2));%dot(p_target, squeeze(min(spec_inf_array))');
        end

        function specinf = specific_information(obj, p_distr, specific_val_dim, specific_val_index)
            % Function to calculate specific information using an probaility distribution
            % This version will be written taking in mind only two variables
            % (p_distr has two dimensions) and later it will expanded to a generic
            % number of dimensions
            % Inputs:
            %   - p_distr: Probability distribution with two dimensions.
            %   - specific_val_dim: Dimension for which specific information is calculated.
            %   - specific_val_index: Index of the specific value in the specified dimension.
            %
            % Output:
            %   - specinf: Specific information in bits.

            d = specific_val_dim;
            v = 1:ndims(p_distr);
            v([1,d]) = v([d,1]);
            bits=1/log(2);
            permuted_prob = permute(p_distr,v); % reshape(permute(p_distr,v),size(p_distr,d),[]);
            marg_prob_specific_dim = sum(permuted_prob, 2:ndims(p_distr));
            marg_prob_specific_val = marg_prob_specific_dim(specific_val_index);
            marg_prob_other_dims = squeeze(sum(permuted_prob, 1));
            specinf = 0;
            for i=1:size(permuted_prob,2)
                if permuted_prob(specific_val_index,i)>0
                    quot = permuted_prob(specific_val_index,i)/marg_prob_specific_val;
                    specinf = specinf +  quot * log(quot/marg_prob_other_dims(i));
                end
            end
            specinf = specinf * bits;
        end

        function atom = calculate_atom(obj, p_distr, target, sources)
            imin = obj.Imin(p_distr, target, sources);
            down_nodes = obj.get_strict_down(sources);
            atom = imin;
            for inodes=1:numel(down_nodes)
                atom = atom - obj.calculate_atom(p_distr, target, down_nodes{inodes});
            end
        end

        function latvals = calculate_latvals(obj, p_distr)
            for node =1:length(obj.lat)
                sources = obj.lat{1,node};
                target = obj.nsources+1;
                obj.pi(node) = obj.calculate_atom(p_distr, target, sources);
            end
        end
        
        function MI = mutualInformation(obj,pxy)
            % Calculate marginal probabilities
            px = sum(pxy, 2); % Marginal probability density of X
            py = sum(pxy, 1); % Marginal probability density of Y

            % Initialize mutual information
            MI = 0;

            % Loop through each value of x and y
            for i = 1:size(pxy, 1)
                for j = 1:size(pxy, 2)
                    if pxy(i,j) ~= 0 && px(i) ~= 0 && py(j) ~= 0
                        MI = MI + pxy(i,j) * log2(pxy(i,j) / (px(i) * py(j)));
                    end
                end
            end
        end
        
        function MIZ = mutualInformationXYZ(obj, pxyz)
            % Calculate marginal probabilities
            pz = sum(sum(pxyz, 1), 2); % Marginal probability density of Z
            pxy = sum(pxyz, 3); % Joint probability density of X and Y

            % Initialize mutual information
            MIZ = 0;

            % Loop through each value of x, y, and z
            for i = 1:size(pxyz, 1)
                for j = 1:size(pxyz, 2)
                    for k = 1:size(pxyz, 3)
                        if pxyz(i,j,k) ~= 0 && pz(k) ~= 0 && pxy(i,j) ~= 0
                            MIZ = MIZ + pxyz(i,j,k) * log2(pxyz(i,j,k) / (pz(k) * pxy(i,j)));
                        end
                    end
                end
            end
        end


        function atom = calculate_atom_mmi(obj, p_distr, target, sources)
            infs = [];
            for s =1:length(sources)
                marg = obj.marginalize(p_distr, [sources, target]);
                infs = [infs obj.mutualInformationLast(marg)];
            end
            
            immi = min(infs);
            down_nodes = obj.get_strict_down(sources);
            atom = immi;
            for inodes=1:numel(down_nodes)
                atom = atom - obj.calculate_atom(p_distr, target, down_nodes{inodes});
            end
        end
        
        function marginal = marginalize(obj, p, keep_dims)
            % p:          N-dimensional probability distribution (array)
            % keep_dims:  vector of dimensions you want to keep (e.g., [1 3])
            
                all_dims = 1:ndims(p);
                sum_dims = setdiff(all_dims, keep_dims);
            
                marginal = p;
                for d = sum_dims
                    marginal = sum(marginal, d);
                end
            
                % Optional: permute to bring kept dims to front
                marginal = permute(marginal, sort(keep_dims));
        end

        function MI = mutualInformationLast(obj, p)
            % Computes I(X;Z) where:
            % - Z is the last dimension of `p`
            % - X is the joint of all other dimensions
            
                % Ensure p is normalized
                p = p / sum(p(:)); 
            
                dims = ndims(p);
                z_dim = dims;
                x_dims = 1:(dims-1);
            
                % Compute marginal p(z)
                pz = sum(p, x_dims);
            
                % Compute marginal p(x) by summing over the last dim
                px = sum(p, z_dim);
            
                % Reshape for broadcasting compatibility
                px_exp = repmat(px, [ones(1, dims-1), size(p, dims)]);
                pz_exp = repmat(reshape(pz, [ones(1, dims-1), numel(pz)]), size(px));
            
                % Avoid division by zero
                valid = p > 0 & px_exp > 0 & pz_exp > 0;
            
                % Compute mutual information
                MI = sum(p(valid) .* log2(p(valid) ./ (px_exp(valid) .* pz_exp(valid))));
            end



    end
end

