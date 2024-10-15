function ii_vector = iiBROJA(pdf_dirty)
%%% Input:
%%%   - pdf_dirty: Probability distribution represented as a multi-dimensional array.
%%%
%%% Output:
%%%   - ii_vector: Integrated information vector containing the minimum self-information value,
%%%                partial information decomposition vector for the permuted distribution, and
%%%                partial information decomposition vector for the original distribution.
if nargin < 1
    msg = 'not enough input arguments.';
    error('iibroja:notEnoughInput', msg);
end

if any(isnan(pdf_dirty), "all")
    msg = 'pdf contains NaNs. Aborting.';
    error('iibroja:NaNInput', msg);
end

if sum(pdf_dirty) == 0
    msg = 'sum of pdf cannot be zero';
    error('iibroja:InvalidInput', msg);
end

if any(pdf_dirty < 0, "all")
    msg = 'negative values in pdf are not allowed';
    error('iibroja:InvalidInput', msg);
end


pid_s = pidBROJA(pdf_dirty);
pdf_choice = permute(pdf_dirty, [3 1 2]);
pid_c = pidBROJA(pdf_choice);
min_si = min(pid_s(1), pid_c(1));

ii_vector = [min_si, pid_c, pid_s];
end