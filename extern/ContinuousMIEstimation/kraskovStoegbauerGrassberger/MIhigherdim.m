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


function miout=MIhigherdim(x,kneig,emb_dim,emb_tau);

% Calculate MI value of any high-dimensional vector - to include also
% dependencies in time use embedding dimension >1 (see[1], Fig 14)
% (rectangular version)
% x....input data mxn   m...channelnummer  n...sampling points  m<<n
% kneig... k nearest neigbor for MI algorithm
% emb_dim... embxding dimension (default is 1, no embedding)
% emb_tau... time-delay, only relevant when emb_dim>1 (default is 1)

%default-values
if ~exist('kneig'), kneig=6; end
if ~exist('emb_dim'), emb_dim=1; end
if ~exist('emb_tau'), emb_tau=1; end

[Nd,N]=size(x);


if Nd>N
    zwsp=x;
    [Nd,N]=size(x);
else
    zwsp=x';
end

save zwspMIhigh.txt zwsp -ASCII


% execute C Programm
[a unout]=unix(['MIhigherdim zwspMIhigh.txt ',num2str(Nd),' ',num2str(emb_dim),' ',num2str(emb_tau),' ',num2str(N),' ',num2str(kneig)]);
miout=str2num(unout);


