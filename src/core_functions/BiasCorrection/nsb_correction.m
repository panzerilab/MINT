function [corrected_v, plugin_v] = nsb_correction(inputs, outputs, corefunc, varargin)
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

new_outputs = outputs;
new_opts = varargin{1};
new_opts.bias = {'plugin'};
for i = 1:length(outputs)
    switch outputs{i}
        case 'I(A;B)'
            new_outputs{i} = 'Insb(A;B)';
    end
end

[corrected_v, plugin_v] = corefunc(inputs, new_outputs, new_opts);

end
