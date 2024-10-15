function cfit_out = cfit_core(p_S4, opts)

pidLattice3 = pid_lattice(3);
pidLattice4 = pid_lattice(4);

node_location3 = {1, 2};
node_location4 = {1, 2, 3};
p_Yt4 = permute(p_S4, [5 2 3 4 1]);
p_S3  = squeeze(sum(p_S4,4));
p_Yt3 = permute(p_S3, [4 2 3 1]);
atoms3 = zeros(2,1);
atoms4 = zeros(2,1);
atoms3(1) =  pidLattice3.calculate_atom(p_S3,  4, node_location3);
atoms3(2) =  pidLattice3.calculate_atom(p_Yt3, 4, node_location3);
           
atoms4(1) =  pidLattice4.calculate_atom(p_S4,  5, node_location4);
atoms4(2) =  pidLattice4.calculate_atom(p_Yt4, 5, node_location4);

cfit_out = min(atoms3) - min(atoms4);
end