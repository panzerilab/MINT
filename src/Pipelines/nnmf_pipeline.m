function W_new = nnmf_pipeline(X, newnDim)
    % Input:
    % X: Input data matrix of size nDim x nTrials (must contain non-negative values)
    % newnDim: The target number of reduced dimensions (newnDim)
    % Output:
    % W_new: Reduced data matrix of size newnDim x nTrials
    
    % Ensure the input matrix X is non-negative
    if any(X(:) < 0)
        error('NNMF requires the input matrix X to contain only non-negative values.');
    end
    
    % Apply NNMF to reduce the dimensionality
    [W, ~] = nnmf(X', newnDim);
    
    % The matrix H contains the reduced representation of X with newnDim dimensions
    W_new = W';  % The new reduced matrix has size newnDim x nTrials
end
