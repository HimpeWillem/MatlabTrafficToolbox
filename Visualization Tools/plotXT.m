function [h1]=plotXT(links,path,load,dt,totT)
%Plots a space time plot
%
%SYNTAX
%   [handle_fig]=plot_speed_xt(startNode,endNode,load,links,timeSteps,ti)
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

%local rename link properties
len=links.length;

coordinatesEnd=cumsum([0;len(path)]');
x=0:length(coordinatesEnd)-1;
x1=0.5:1:length(coordinatesEnd)-1;
coordinates = interp1(x,coordinatesEnd,x1)';
plotLoad=load(path,:);

handle_fig = figure('Units','pixels');
hold on
surf(0:dt:totT*dt,coordinates,plotLoad,'facecolor','texturemap','LineStyle','none')
view(2)
axis('tight') 
xlabel('Time [hr]','FontSize',12)
ylabel('Distance [km]','FontSize',12)
caxis([0 max(max(plotLoad))]);
colorbar('EastOutside');
my_jet = jet(512);
my_jet(1:250,:)=[];
colormap(my_jet);
grid on
title('Space - Time graph','FontSize',14,'fontweight','b')
% set(handle_fig, 'Position', [100, 100, 400, 300]);
