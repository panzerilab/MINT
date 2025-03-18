
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
% Check Outputslist

function bias = pt(C_ptr, Rtot, N)
    % Function to compute the bias as described in Panzeri & Treves (1996)
    % Input:
    %   C_ptr - Array of values
    %   Rtot - Total number of bins
    %   N - Sample size
    % Output:
    %   bias - Computed bias value

    % Number of non-zero C values
    NnonZeroCvalues = sum(C_ptr(:) > 0);
    
    % Initialize PnonZero array
    PnonZero = zeros(1, NnonZeroCvalues);
    
    % Read non-zero probability values and compute Rplugin
    index = 1;
    for i = 1:Rtot
        if C_ptr(i) > 0
            PnonZero(index) = C_ptr(i); % / N;
            index = index + 1;
        end
    end
    
    % Rplugin is the number of non-null probability values
    Rplugin = NnonZeroCvalues;
    
    if Rplugin < Rtot
        % Initial value for Rexpected
        Rexpected = Rplugin;
        for i = 1:NnonZeroCvalues
            Rexpected = Rexpected - (1 - PnonZero(i))^N;
        end
        
        % Initial values for deltaRprevious and deltaR
        deltaRprevious = Rtot;
        deltaR = abs(Rplugin - Rexpected);
        
        R = Rplugin;
        while deltaR < deltaRprevious && R < Rtot
            R = R + 1;
            
            gamma = (R - Rplugin) * (1 - (N / (N + Rplugin))^(1 / N));
            
            Rexpected = 0;
            % Occupied bins
            for i = 1:NnonZeroCvalues
                Pbayes = (PnonZero(i) * N + 1) / (N + Rplugin) * (1 - gamma);
                Rexpected = Rexpected + 1 - (1 - Pbayes)^N;
            end
            
            % Non-occupied bins
            Pbayes = gamma / (R - Rplugin);
            Rexpected = Rexpected + (R - Rplugin) * (1 - (1 - Pbayes)^N);
            
            deltaRprevious = deltaR;
            deltaR = abs(Rplugin - Rexpected);
        end
        
        R = R - 1;
        if deltaR < deltaRprevious
            R = R + 1;
        end
    else
        R = Rtot;
    end
    
    % Estimating bias
    bias = -(R - 1) / (2 * N * log(2));
end
