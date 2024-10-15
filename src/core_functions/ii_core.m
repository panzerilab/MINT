function atoms = ii_core(p_src, opts)
p_crs = permute(p_src, [3 2 1]);
atoms = zeros(2,1);
atoms_1 = feval(opts.function, p_src);
atoms_2 = feval(opts.function, p_crs);
atoms(1) = atoms_1(1);
atoms(2) = atoms_2(1);
end
