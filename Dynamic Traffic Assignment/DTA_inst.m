function [cvn_up,cvn_down,TF] = DTA_inst(nodes,links,origins,destinations,ODmatrix,dt,totT,rc_dt)

%Method of successive averages for calculating deterministic user
%equilibrium
%
%SYNTAX
%   [flows] = MSA_exercise(odmatrix,nodes,links)
%
%DESCRIPTION
%   returns the flow on each link in the deterministic user equilibrium as
%   calculated by the method of successive averages
%
%INPUTS
%   odmatrix: static origin/destination matrix
%   nodes: table with all the nodes in the network.
%   links: table with all the links in the network
    
% This file is part of the matlab package for dynamic traffic assignments 
% developed by the KULeuven. 
%
% Copyright (C) 2016  Himpe Willem, Leuven, Belgium
%
% This program is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% any later version.
% 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
% More information at: http://www.mech.kuleuven.be/en/cib/traffic/downloads
% or contact: willem.himpe {@} kuleuven.be

%initilization
totNodes = size(nodes.id,1);
totLinks = size(links.toNode,1);
totDest = length(destinations);
%cumulative vehicle numbers (cvn) are stored on both upstream and
%dowsntream link end of each link for every time slice
cvn_up = zeros(totLinks,totT+1,totDest);
cvn_down = zeros(totLinks,totT+1,totDest);
    
%calculate new flows
[cvn_up,cvn_down,TF] = LTM_MC_INST(nodes,links,origins,destinations,ODmatrix,dt,totT,rc_dt);

%update costs
[simTT] = cvn2tt(sum(cvn_up,3),sum(cvn_down,3),dt,totT,links);

%Compute the convergence gap
%using the Null option no turning fractions are calculated
[~,gap] = allOrNothingTF(nodes,links,destinations,simTT,cvn_up,dt,totT,rc_dt,'Null');

%Check for number of iterations until convergence
disp(['Instantaneous route choice: No Iterations', ' Gap: ', num2str(gap)]);

end


