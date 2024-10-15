function atoms = fit_core(p_S, opts)
node_location = {1,2};
pidLattice = pid_lattice(3);
p_Yt = permute(p_S, [4 2 3 1]);
atoms = zeros(2,1);
atoms(1) = pidLattice.calculate_atom(p_S,  4, node_location);
atoms(2) = pidLattice.calculate_atom(p_Yt, 4, node_location);
end