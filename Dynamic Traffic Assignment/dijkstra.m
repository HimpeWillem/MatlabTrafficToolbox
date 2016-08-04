function [parents, distance, path] = dijkstra(matrix,source,varargin)

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

%initialize
n = size(matrix, 2);
parents = zeros(n,1);
distance= inf(n,1);
visited = false(n,1);

%use priority queue for fast ordering and lookup of minimal distance
comp = javaObject('PairComparator');
nodesQueue = java.util.PriorityQueue(max(10,ceil(n/10)),comp);

%add the source node the list
nodesQueue.add([0, source]);

if nargin==2
    %if no destination node is provided the algorithm runs until all nodes 
    %are labeled
    maxnum = numel(nodesQueue.peek());
    
    while numel(nodesQueue.peek()) > 0 %nodesQueue.size()>0
        %remove closest node
        elem = nodesQueue.poll();
        u = elem(2);%node
        k = elem(1);%travel cost
        %if it is already visited ignore it
        if visited(u)
            continue;
        end
        %fix the distance
        distance(u)=k;
        visited(u)=true;
        %go over all unvisited neighbours
        [~,nv,s] = find(matrix(u,:));
        for i = 1:length(nv)
            v=nv(i);
            if~visited(v)
                temp = s(i) + distance(u);
                if temp < distance(v)
                    %if a new value is found it is added to the list
                    nodesQueue.add([temp, v]);
                    distance(v) = temp;
                    parents(v) = u;
                end
            end
        end
    end
    path=[];
else
    %if a destination node is provided the algorithm runs until that node 
    %is labeled
    destination = varargin{1};
    u=source;  
    while u~=destination && numel(nodesQueue.peek()) > 0
        %remove closest node
        elem = nodesQueue.poll();
        u = elem(2);%node
        k = elem(1);%travel cost
        %if it is already visited ignore it
        if visited(u)
            continue;
        end
        %fix the distance
        distance(u)=k;
        visited(u)=true;
        %go over all unvisited neighbours
        [~,nv,s] = find(matrix(u,:));
        for i = 1:length(nv)
            v=nv(i);
            if~visited(v)
                temp = s(i) + distance(u);
                if temp < distance(v)
                    %if a new value is found it is added to the list
                    nodesQueue.add([temp, v]);
                    distance(v) = temp;
                    parents(v) = u;
                end
            end
        end
    end
        
    %get the nodes that forms the shortest path
    path = dijkstra_getpath(source,destination,parents);
end

    
    
    