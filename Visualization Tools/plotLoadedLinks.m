function [handle_fig,handle_rect,handle_txt] = plotLoadedLinks(nodes,links,load,show_labels,fig_num,scale,maxLoad)
%Plots the link loads on a network
%
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
%
%
%SYNTAX
%   [handle_fig,handle_rect,handle_txt] = plotLoadedLinks(nodes,links,load,show_labels,fig_num,scale,maxLoad)
%
%DESCRIPTION
%   plots loads on link in the network and returns figure handles of the loaded network
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
strN = links.fromNode;
endN = links.toNode;
x=nodes.xco;
y=nodes.yco;

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
    
    %Setup figure
    margX=0.1*(max(x)-min(x))+100*eps;
    margY=0.1*(max(y)-min(y))+100*eps;
    axis([min(x)-margX max(x)+margX min(y)-margY max(y)+margY]);
    colorbar('EastOutside');

else
    handle_fig = figure(fig_num);
    hold on;
end

if isempty(maxLoad)
    maxc=max(load);
else
    maxc=maxLoad;
end

% handle_ax=axes;

%make a rectangular object for each link
upX=x(strN);
downX=x(endN);
upY=y(strN);
downY=y(endN);

%set scale
if isempty(scale)
    scale = 1/maxc*0.25*sqrt((max(x)-min(x))^2+(max(y)-min(y))^2)/length(upX)^(1/2);
end

vx=downX-upX;
vy=downY-upY;
vl=sqrt(vx.^2+vy.^2);
vx=vx./vl;
vy=vy./vl;

sc=scale;
sc=eps+scale*load';

xrec=[upX';upX'+sc.*vy';downX'+sc.*vy';downX'];
yrec=[upY';upY'-sc.*vx';downY'-sc.*vx';downY'];


%set the colours
ctemp=hsv(128);
cmap=colormap(ctemp(50:-1:1,:));
% cmap=colormap(handle_ax,ctemp(50:-1:1,:));


minc=0;%possible one could also use the minimal positive value of the load  %max(0,min(load));

crec=cmap(ceil(49*(load'-minc+eps)/(maxc-minc+eps))',:);
caxis([minc maxc]);

%visualize all loads
handle_rect=patch(xrec,yrec,load');
set(handle_rect,'FaceColor','flat','FaceVertexCData',crec);
warning off verbose
% colorbar('EastOutside');

handle_txt=[];
if show_labels
    cx=[x(strN)+x(endN)]/2;
    cy=[y(strN)+y(endN)]/2;
    id=[1:length(cx)]';
    bl=(upX >= downX & upY <= downY);
    tl=(upX < downX & upY < downY);
    br=(upX >= downX & upY > downY);
    tr=(upX < downX & upY >= downY);
    t1=text(cx(bl),cy(bl),num2str(load(bl)),'Color',[0 0 0],'VerticalAlignment','Bottom','HorizontalAlignment','Left','FontWeight','bold','Clipping','on','hittest','off');
    t2=text(cx(tl),cy(tl),num2str(load(tl)),'Color',[0 0 0],'VerticalAlignment','Top','HorizontalAlignment','Left','FontWeight','bold','Clipping','on','hittest','off');
    t3=text(cx(br),cy(br),num2str(load(br)),'Color',[0 0 0],'VerticalAlignment','Bottom','HorizontalAlignment','Right','FontWeight','bold','Clipping','on','hittest','off');
    t4=text(cx(tr),cy(tr),num2str(load(tr)),'Color',[0 0 0],'VerticalAlignment','Top','HorizontalAlignment','Right','FontWeight','bold','Clipping','on','hittest','off');       
	handle_txt=[t1;t2;t3;t4];
end

% set(handle_fig, 'Position', [100, 100, 500, 300]);
hold off;

end