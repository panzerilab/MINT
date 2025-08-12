function p_ts = create_prob_ts(p_distr, dims)
    % dims = [target, source] where target must be scalar
    % Output is 2D: (|target|) x (|sources_combined|)
    
    % Step 1: marginalize over all other dimensions
    keep_dims = dims;
    all_dims = 1:ndims(p_distr);
    marginalize_dims = setdiff(all_dims, keep_dims);
    p_ts = p_distr;
    
    for d = sort(marginalize_dims, 'descend')
        p_ts = sum(p_ts, d);
    end

    % Step 2: permute so target is first
    [~, perm] = ismember(dims, all_dims);
    
    % Fix: ensure perm has correct length
    if numel(perm) < ndims(p_ts)
        perm = [perm, setdiff(1:ndims(p_ts), perm, 'stable')];
    end
    
    p_ts = permute(p_ts, perm);

    % Step 3: reshape into 2D: target x sources_combined
    target_size = size(p_ts, 1);
    source_size = numel(p_ts) / target_size;
    p_ts = reshape(p_ts, target_size, source_size);
end
