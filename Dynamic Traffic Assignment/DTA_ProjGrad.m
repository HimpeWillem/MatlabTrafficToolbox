function [cvn_up,cvn_down,TF] = DTA_ProjGrad(nodes,links,origins,destinations,ODmatrix,dt,totT,rc_dt,maxIt,alpha)

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

%setup the output figure
h = figure;
semilogy(0,NaN);
start_time = cputime;

%Maximum number of iterations
if isempty(maxIt)
    maxIt = 20;
end

%initilization
totNodes = size(nodes.id,1);
totLinks = size(links.toNode,1);
totDest = length(destinations);
%cumulative vehicle numbers (cvn) are stored on both upstream and
%dowsntream link end of each link for every time slice
cvn_up = zeros(totLinks,totT+1,totDest);
cvn_down = zeros(totLinks,totT+1,totDest);

%Initialize the iteration numbering
it = 0;

%Initialize the aggragation
rc_agg = 'last'; %[first/middle/last]

%initialize the travel cost
[simTT] = cvn2tt(sum(cvn_up,3),sum(cvn_down,3),dt,totT,links);

%initialize the gap function
gap_dt = inf;
TF = [];

%initialize the turning fractions (all options [first/middle/last] give the same result)
TF_new = allOrNothingTF(nodes,links,destinations,simTT,cvn_up,dt,totT,rc_dt,rc_agg);

%MAIN LOOP: iterate until convergence is reached or maximum number of
%iterations is reached
while it < maxIt && gap_dt > 10^-60
    it = it+1;
    
    %calculate the update step
    if isempty(TF)
        TF = TF_new;        
    else
        for n = 1:totNodes
            for t = 1:totT
                for d= 1:totDest
                    TF{n,t,d} = TF{n,t,d} + updateTF{n,t,d};
                end
            end
        end
    end
   
    %calculate new flows
    [cvn_up,cvn_down] = LTM_MC(nodes,links,origins,destinations,ODmatrix,dt,totT,TF);
    
    %update costs
    [simTT] = cvn2tt(sum(cvn_up,3),sum(cvn_down,3),dt,totT,links);
    
    %Compute new turning fractions via all-or-nothing assignment
    %and compute convergence gap
    [updateTF,gap_dt,gap_rc] = projGradTF(nodes,links,destinations,simTT,cvn_up,dt,totT,rc_dt,rc_agg,TF,alpha);
    
    %plot convergence in function of computation time
    figure(h)
    hold on
    time=cputime-start_time;
    a=semilogy(time,gap_dt,'r.');
    b=semilogy(time,gap_rc,'ob');
    legend([a,b],'gap based on simulation interval','gap based on route choice interval');
    
    sp=[TF{2,:,1}];
    figure(10);plot(dt*[0:totT-1],sp(1:2:end),'r',dt*[0:totT-1],sp(2:2:end),'b');
    title(['iteration: ',num2str(it)]); 
    pause(0.01);
end
%Check for number of iterations until convergence
if it >= maxIt
    disp(['Maximum Iteration limit reached: ', num2str(maxIt), ' Gap: ', num2str(gap_dt)]);
else
    disp(['Convergence reached in iteration ', num2str(it), ' Gap: ', num2str(gap_dt)]);
end

end


