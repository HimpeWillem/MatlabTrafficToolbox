function [handle_fig,handle_sctr,handle_txt] = plotLoadedNodes(nodes,links,load,show_label,fig_num,scale,maxLoad)
%Plots the node loads on a network
%
%SYNTAX
%   [handle_fig,handle_scatter,handle_txt] = plotLoadedNodes(nodes,links,load,show_label,fig,scale,maxLoad)
%
%DESCRIPTION
%    plots loads on link in the network and returns figure handles of the loaded network
%
%INPUTS
%   nodes: list of all the nodes in the network.
%   Each entry of the list represents one node. Each node is a structure that
%   has at least a node ID and an x and y coordinate of the node
%   links: list of all the links in the network
%   Each entry of the list represents one link. Each link is a structure that
%   has at least a link ID and an upstream and downstream node.
%   load: list of the load on each link
%   show_text: flag for visualizing labels
%   fig: parameter for the figure handle
%   maxLoad: scale parameter used for determining the loads

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
x=nodes.xco;
y=nodes.yco;
strN = links.fromNode;
endN = links.toNode;

totN = length(x);

%open figure
if isempty(fig_num)
    handle_fig = figure;
    hold on;
    plot(x,y,'.','MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0 0 0]);
    
    x_temp = zeros(length(strN)*3,1);
    x_temp(1:3:end) = x(strN);
    x_temp(2:3:end) = x(endN);
    x_temp(3:3:end) = NaN;
    
    y_temp = zeros(length(strN)*3,1);
    y_temp(1:3:end) = y(strN);
    y_temp(2:3:end) = y(endN);
    y_temp(3:3:end) = NaN;
    
    plot(x_temp, y_temp,'Color',[0 0 0]); %[0.9 0.7 0]
else
    handle_fig = figure(fig_num);
    hold on;
end
% handle_ax=axes;

%set scale
if isempty(scale)
    scale = 1/max(load)*250*sqrt((max(x)-min(x))^2+(max(y)-min(y))^2)/length(strN)^1.5;
end

%make a point object for each link
sc=scale*load';


%set the colours
ctemp=hsv(128);
cmap=colormap(ctemp(128:-1:78,:));
% cmap=colormap(handle_ax,ctemp(128:-1:78,:));

minc=0;%possible one could also use the minimal positive value of the load  %max(0,min(load));

if isempty(maxLoad)
    maxc=max(load);
else
    maxc=maxLoad;
end

crec=cmap(ceil(49*(load'-minc+10*eps)/(maxc-minc+10*eps))',:);

%plot all loads
handle_sctr=scatter(x,y,eps+sc,crec,'filled');

caxis([minc maxc]);
colorbar('EastOutside');
% colorbar(handle_ax,'WestOutside');

handle_txt=[];
if show_label
    
    a = load>=0;
    handle_txt=text(x(a),y(a),num2str(load(a)),'Color',[0 0 0],'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontWeight','bold');
end

%Setup figure
margX=0.1*(max(x)-min(x))+100*eps;
margY=0.1*(max(y)-min(y))+100*eps;
axis([min(x)-margX max(x)+margX min(y)-margY max(y)+margY])


% set(handle_fig, 'Position', [100, 100, 500, 300]);
hold off;

end