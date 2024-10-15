function [corrected_v, naive_v] = shufflecorrection(inputs, outputs, corefunc, varargin)

new_outputs = outputs;
new_opts = varargin{1};
new_opts.bias = {'naive'};
for i = 1:length(outputs)
    switch outputs{i}
        case 'I(A;B)'
            new_outputs{i} = 'Ish(A;B)';
        case 'coI(A;B)'
            new_outputs{i} = 'coIsh(A;B)';
        case 'Ic(A;B)'
            new_outputs{i} = 'Icsh(A;B)';
        case 'Icd(A;B)'
            new_outputs{i} = 'Icdsh(A;B)';
    end
end

[corrected_v, naive_v] = corefunc(inputs, new_outputs, new_opts);

end
