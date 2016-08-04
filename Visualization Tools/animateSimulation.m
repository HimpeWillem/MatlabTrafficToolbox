function animateSimulation(nodes,links,values,timeSteps,frate)
%Animate a simulation
%
%
%
%SYNTAX
%   animateSimulation(nodes,links,values,timeSteps,option)
%
%DESCRIPTION
%   simulate the loads on the network
%
%INPUTS
%   nodes: list of all the nodes in the network.
%   Each entry of the list represents one node. Each node is a structure that
%   has at least a node ID and an x and y coordinate of the node
%   links: list of all the links in the network
%   Each entry of the list represents one link. Each link is a structure that
%   has at least a link ID and an upstream and downstream node.
%   values: time dependent loads
%   timeSteps: simulation time of each step
%   scale: scale
%   option:

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

display('press space bar to continue');
[handle_fig,handle_txt]=plotNetwork(nodes,links,false,11);

%set axis
strN = links.fromNode;
endN = links.toNode;
load = max(values,[],2);
x=nodes.xco;
y=nodes.yco;
upX=x(strN);
downX=x(endN);
upY=y(strN);
downY=y(endN);
 maxc=max(max(values));
scale = 1/maxc*0.25*sqrt((max(x)-min(x))^2+(max(y)-min(y))^2)/length(upX)^(1/2);
vx=downX-upX;
vy=downY-upY;
vl=sqrt(vx.^2+vy.^2);
vx=vx./vl;
vy=vy./vl;
sc=scale;
sc=eps+scale*load';
xrec=[upX';upX'+sc.*vy';downX'+sc.*vy';downX'];
yrec=[upY';upY'-sc.*vx';downY'-sc.*vx';downY'];
margX=0.1*(max(x)-min(x))+100*eps;
margY=0.1*(max(y)-min(y))+100*eps;
axis([min(min(xrec))-margX max(max(xrec))+margX min(min(yrec))-margY max(max(yrec))+margY])

pause()
delete(handle_txt);
colorbar('EastOutside');
for t=1:size(timeSteps,2)
    title(num2str(timeSteps(t)));
    [handle_fig,handle_rect,handle_txt] = plotLoadedLinks(nodes,links,values(:,t),false,handle_fig,[],max(max(values(:,:))));

    if isinf(frate)
        pause();
    else
        pause(1/frate)
    end
    if ~isempty(handle_rect)
        try
            delete(handle_rect);
        end
    end
    if~isempty(handle_txt)
        try
            delete(handle_txt);
        end
    end
    if ~ishandle(handle_fig)
       return;
    end
end
   
end