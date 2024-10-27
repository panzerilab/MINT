%%% Description:
%%% This class implements a redundancy lattice, a mathematical structure used to model redundancy relationships between sets of variables. It is particularly designed for use in the context of probabilistic information theory.
%%% 
%%% Properties:
%%% - lat: A cell array representing the lattice nodes.
%%% - red: Redundancy values associated with lattice nodes.
%%% - pi: Probability distribution array.
%%% - pdf: Probability density function associated with the lattice nodes.
%%% 
%%% Methods:
%%% Constructor: pid_lattice(nsources)
%%%    - Initializes the pid_lattice object with an empty lattice and generates lattice nodes based on the given number of variables (nsources).
%%% 
%%% power_set(obj, var_list)
%%%    - Computes the power set of a given set of variables using bitwise operations in MATLAB.
%%% 
%%% is_node_red(obj, x)
%%%    - Checks if a given collection of nodes x represents a valid node in the redundancy lattice.
%%% 
%%% node_issubsetany(obj, x, i, j)
%%%    - Checks if node x[i] is a subset of any other node x[j] within a given collection.
%%% 
%%% node_issubset(obj, x, y)
%%%    - Determines if a node x precedes another node y in the redundancy lattice.
%%% 
%%% node_issame(obj, x, y)
%%%    - Checks if two nodes x and y are identical in terms of their elements.
%%% 
%%% get_down(obj, node)
%%%    - Returns the set of lattice nodes that are below the given node in the lattice.
%%% 
%%% get_strict_down(obj, node)
%%%    - Returns the set of lattice nodes that are strictly below the given node in the lattice.
%%% 
%%% get_red(obj, node)
%%%    - Computes the redundancy value associated with the given node.
%%% 
%%% Imin(obj, p_distr, target, sources)
%%%     - Computes the minimum specific information between a target variable and multiple source variables, given a probability distribution.
%%% 
%%% specific_information(obj, p_distr, specific_val_dim, specific_val_index)
%%%     - Calculates specific information between two variables based on a probability distribution.
%%% 
%%% Note: The class employs various methods to handle lattice nodes, check redundancy relationships, and calculate information-theoretic measures.
