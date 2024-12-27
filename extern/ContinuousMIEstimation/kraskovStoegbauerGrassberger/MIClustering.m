% Copyright 2009 Alexander Kraskov, Harald Stoegbauer, Ralph G. Andrzejak, Peter Grassberger
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
% A. Kraskov, H. Stogbauer, R. G. Andrzejak, and P. Grassberger, 
% Hierarchical clustering using mutual information.
% EUROPHYSICS LETTERS 70 (2): 278-284,  2005
%
% in your published research.


function MIClustering(x,kneig,algo)

% hierarchical clustering using the grouping property of mutual information
% output: dendrogram 
% x....input data mxn   m...channelnummer  n...sampling points  m<<n
% kneig... k nearest neigbor for MI algorithm
% algo ... 1=cubic  2=rectangular

%default-values
if ~exist('kneig'), kneig=6; end
if ~exist('algo'), algo=2; end

% check input data if format is correct
[Nd,N]=size(x);
if Nd>N
    x=x';
    [Nd,N]=size(x);
end

zwsp=x';
save zwspMIC.txt zwsp -ASCII

% execute C Programm
[a b]=unix(['MIClustering zwspMIC.txt ',num2str(Nd),' ',num2str(N),' ',num2str(kneig),' ',num2str(algo)]);
clutr=str2num(b);

% Plot Output
figure;
dendrogram(clutr(:,[1,2,5]));
ylabel('mutual information');
title('Hierarchical clustering using the grouping property of mutual information');


