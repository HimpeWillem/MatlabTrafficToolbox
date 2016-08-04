function [TF_P,arr_map] = TF_init(node_prop,links,destinations,dt,totT)
%#codegen

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

%based on free flow travel times

tt_base=links.length./links.freeSpeed;
totDest = length(destinations);
totNod = length(node_prop.nbIncomingLinks);
totTF = length(node_prop.outgoingLinksTF);

TF_P = zeros(totTF,totDest,totT);
distance_destination = zeros(totNod,1);
temp_distance = zeros(totNod,1);
parent = zeros(totNod,1);

arr_map = zeros(totNod,totT+1,totDest);

%simple proportion matrix
%in free flow conditions (this is used in the first DTA iteration)
for d_index=1:totDest
    tt=tt_base;
    for n_index=1:totDest
        if n_index~=d_index
            n = destinations(n_index);
            tt(node_prop.incomingLinksList(node_prop.positionFirstIn(n):(node_prop.positionFirstIn(n)+node_prop.nbIncomingLinks(n)-1)))=inf;
        end
    end
    
    d=destinations(d_index);
    %dijkstra tree
    distance_destination(:,1) = inf;    % it stores the shortest distance between each node and the source node;
    temp_distance(:,1) = inf;
    TF_P(:,d_index,:)=0;
    parent(:,1) = 0;
    distance_destination(d) = 0;
    temp_distance(d) = 0;
    
    for i = 1:totNod
        [time, u] = min(temp_distance(:,1));       % it starts from node with the shortest distance to the source;
        temp_distance(u) = inf;                    % mark it as visited;
        
        incomingLinks = node_prop.incomingLinksList(node_prop.positionFirstIn(u):(node_prop.positionFirstIn(u)+node_prop.nbIncomingLinks(u)-1));
        
        for l_index=1:node_prop.nbIncomingLinks(u)            % for each neighbors of node u;
            l=incomingLinks(l_index);
            s=links.fromNode(l);
            if (tt(l) + time < distance_destination(s))
                distance_destination(s) = time + tt(l);   % update the shortest distance when a shorter path is found;
                temp_distance(s) = distance_destination(s);
                parent(s) = u;                                     % update its parent;
            end;
        end;
    end
    
    lout=find(parent(links.fromNode)==links.toNode);
    for l=1:length(lout)
        TF_P(node_prop.outgoingLinksTF==lout(l),d_index,:)=1;
    end
    
%     out_link_P(parent(links.startNodes)==links.endNodes,d_index,:)=1;
    arr_map(:,:,d_index)=distance_destination*ones(1,totT+1)+dt*ones(totNod,1)*[0:1:totT];
end