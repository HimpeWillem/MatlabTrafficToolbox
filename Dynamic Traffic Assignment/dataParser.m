function [links,node_prop] = dataParser(links,nodes,origins,destinations,dt)

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

%Setup nodes
node_prop=[];

%connectivity (incoming and outgoing links)
totNodes = length(nodes.id);

for n=1:totNodes
    node(n).incomingLinks = find([links.toNode]==n)';
    node(n).outgoingLinks = find([links.fromNode]==n)';
end

%connectivity of nodes with incoming links and outgoing links
node_prop.incomingLinksList=[node.incomingLinks]';
node_prop.outgoingLinksList=[node.outgoingLinks]';

node_prop.nbIncomingLinks=zeros(totNodes,1);
node_prop.nbOutgoingLinks=zeros(totNodes,1);

node_prop.positionFirstIn=zeros(totNodes,1);
node_prop.positionFirstOut=zeros(totNodes,1);
node_prop.positionFirstIn(1)=1;
node_prop.positionFirstOut(1)=1;
node_prop.nbIncomingLinks(1)=length(node(1).incomingLinks);
node_prop.nbOutgoingLinks(1)=length(node(1).outgoingLinks);


for n=2:totNodes
    node_prop.nbIncomingLinks(n)=length(node(n).incomingLinks);
    node_prop.nbOutgoingLinks(n)=length(node(n).outgoingLinks);
    node_prop.positionFirstIn(n)=node_prop.positionFirstIn(n-1)+node_prop.nbIncomingLinks(n-1);
    node_prop.positionFirstOut(n)=node_prop.positionFirstOut(n-1)+node_prop.nbOutgoingLinks(n-1);
end

%Turning Fractions setup
node_prop.nbTF = (max(1,node_prop.nbIncomingLinks)).*(max(1,node_prop.nbOutgoingLinks));

node_prop.incomingLinksTF = zeros(sum(node_prop.nbTF),1);
node_prop.incomingLinksTFindex = zeros(sum(node_prop.nbTF),1);
node_prop.outgoingLinksTF = zeros(sum(node_prop.nbTF),1);
node_prop.outgoingLinksTFindex = zeros(sum(node_prop.nbTF),1);

node_prop.positionFirstTF = zeros(totNodes,1);

lpos=1;
for n=1:totNodes
    if any(n==origins)
        node_prop.positionFirstTF(n)=lpos;
        for lout=1:max(1,node_prop.nbOutgoingLinks(n))
            node_prop.incomingLinksTF(lpos) = 0;
            if node_prop.nbOutgoingLinks(n)>0
                node_prop.outgoingLinksTF(lpos) = node_prop.outgoingLinksList(node_prop.positionFirstOut(n)+lout-1);
                node_prop.outgoingLinksTFindex(lpos) = lout;
            else
                node_prop.outgoingLinksTF(lpos) = 0;
            end
            lpos=lpos+1;
        end
    elseif node_prop.nbIncomingLinks(n)>0 || node_prop.nbOutgoingLinks(n)>0
        node_prop.positionFirstTF(n)=lpos;
        for lin=1:max(1,node_prop.nbIncomingLinks(n))
            for lout=1:max(1,node_prop.nbOutgoingLinks(n))
                if node_prop.nbIncomingLinks(n)>0
                    node_prop.incomingLinksTF(lpos) = node_prop.incomingLinksList(node_prop.positionFirstIn(n)+lin-1);
                    node_prop.incomingLinksTFindex(lpos) = lin;
                else
                    node_prop.incomingLinksTF(lpos) = 0;
                end
                if node_prop.nbOutgoingLinks(n)>0
                    node_prop.outgoingLinksTF(lpos) = node_prop.outgoingLinksList(node_prop.positionFirstOut(n)+lout-1);
                    node_prop.outgoingLinksTFindex(lpos) = lout;
                else
                    node_prop.outgoingLinksTF(lpos) = 0;
                end
                lpos=lpos+1;
            end
        end
    end
end

% Setup links

%spillbackspeed
links.ws = - links.capacity./(links.kJam - links.capacity./links.freeSpeed);

%initialization for ILTM
links.vf_index = floor((links.length./links.freeSpeed)/dt);
links.vf_ratio = links.vf_index-(links.length./links.freeSpeed)/dt+1;
links.vf_index = -links.vf_index-1;
links.vw_index = floor((-links.length./links.ws)/dt);
links.vw_ratio = links.vw_index-(-links.length./links.ws)/dt+1;
links.vw_index = -links.vw_index-1;

%links towards destinations should have infinit storage capacity (because
%they are sinks) if one tries to model a restricted source is should be
%modelled as a constraint on the inflow capacity of the link (not the jam
%densitity)
for d_index=1:length(destinations)
    d=destinations(d_index);
    for l_index=1:node_prop.nbIncomingLinks(d)
        l=node_prop.incomingLinksList(node_prop.positionFirstIn(d)+l_index-1);
        links.kJam(l)=inf;
    end
end

end