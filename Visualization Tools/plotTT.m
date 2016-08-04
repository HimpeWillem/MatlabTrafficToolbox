function [handle_fig,instTT,forwTT,backTT]=plotTT(links,path,TT,dt,totT)
%Plots a travel time plot
%
%
%
%SYNTAX
%   [handle_fig]=plotTT(links,path,TT,dt,totT)
%
%DESCRIPTION
%   simulate the loads on the network
%
%INPUTS
%   startNode: start node
%   endNode: end node
%   load: the load 
%   nodes: list of all the nodes in the network.
%   Each entry of the list represents one node. Each node is a structure that
%   has at least a node ID and an x and y coordinate of the node
%   links: list of all the links in the network
%   Each entry of the list represents one link. Each link is a structure that
%   has at least a link ID and an upstream and downstream node.
%   timesteps: time steps of the simulation
%   ti: title

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


handle_fig = figure('Units','pixels');
xlabel('Time [hr]','FontSize',12);
ylabel('Travel Time [hr]','FontSize',12);
title('Travel time graph','FontSize',14,'fontweight','b');
hold on;

timeSteps=0:dt:totT*dt;
instTT=sum(TT(path,:),1);
plot(timeSteps,instTT,'g')

forwTT=timeSteps;
backTT=timeSteps;
for l=path
    forwTT=forwTT+interp1(timeSteps,TT(l,:),forwTT);
end
forwTT=forwTT-timeSteps;
forwTT(isnan(forwTT))=min(instTT);
plot(timeSteps,forwTT,'b');

for l=path(end:-1:1)
    [x,un]=unique(timeSteps+TT(l,:),'last');
    arTT = interp1(x,TT(l,un),timeSteps);
    arTT(isnan(arTT)) = TT(l,1);
    backTT=backTT-interp1(timeSteps,arTT,backTT);
end
backTT=timeSteps-backTT;
backTT(isnan(backTT))=min(instTT);
plot(timeSteps,backTT,'r');

plot(timeSteps,sum(links.length(path)./links.freeSpeed(path))*ones(1,totT+1),'--k')
grid on
legend('instantaneous','at departure','at arrival (experienced)','free flow travel time') 
axis([timeSteps(1),timeSteps(end),max(0,min(sum(TT(path,:),1)-0.025)),max(sum(TT(path,:),1))+0.025]);
return
