% Copyright 2009 Alexander Kraskov, Harald Stoegbauer, Peter Grassberger
%-----------------------------------------------------------------------------------------
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
% You should receive a copy of the GNU General Public License
% along with this program.  See also <http://www.gnu.org/licenses/>.
%----------------------------------------------------------------------------------------- 
% Contacts:
%
% Harald Stoegbauer <h.stoegbauer@gmail.com>
% Alexander Kraskov <alexander.kraskov@gmail.com>
%-----------------------------------------------------------------------------------------
% Please reference
% 
% A. Kraskov, H. Stogbauer, and P. Grassberger,
% Estimating mutual information.
% Phys. Rev. E 69 (6) 066138, 2004
%
% in your published research.


function miout=MIxnyn_matlab(x,y,kneig, pathToSave);

% Calculate MI value between 2 vector of any dimension (rectangular
% version)
% x....input data mxn   m...channelnummer  n...sampling points  m<<n
% kneig... k nearest neigbor for MI algorithm


%default-values
if ~exist('kneig'), kneig=6; end


% check input data if format is correct
[Ndx,Nx]=size(x);
if Ndx>Nx
    x=x';
    [Ndx,Nx]=size(x);
end
[Ndy,Ny]=size(y);
if Ndy>Ny
    y=y';
    [Ndy,Ny]=size(y);
end

if Nx~=Ny
    if Nx>Ny
        N=Ny;
    else
        N=Nx;
    end
    fprintf('The two input vectors must have the same length !!!!');
    fprintf('Caluculation using the %d datapoints',N);
    
else
    N=Nx;    
end


% save data for C-Programm
zwsp=[x;y]';
save zwspMIxnyn.txt zwsp -ASCII


% execute C Programm
path = [pathToSave '/zwspMIxnyn.txt'];
[miout] = MIxnyn(path,Ndx, Ndy, N,kneig);


