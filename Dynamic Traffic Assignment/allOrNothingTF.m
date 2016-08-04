function [TF,gap_dt,gap_rc] = allOrNothingTF(nodes,links,destinations,simTT,cvn_up,dt,totT,rc_dt,rc_agg)

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

totDest = length(destinations);
totNodes = length(nodes.id);
totLinks = length(links.id);
strN = links.fromNode;
endN = links.toNode;

if isempty(cvn_up)
    cvn_up=zeros(totLinks,totT+1,totDest);
end

TF = num2cell(ones(size(nodes.id,1),totT,totDest));
timeSteps = dt*[0:1:totT];
timeRC = rc_dt*[0:1:totT];
timeRC(timeRC>timeSteps(end))=[];

gap = zeros(totLinks,totT+1);
gap_dt = 0;
gap_rc = 0;
act_t = false(1,totT+1);
gVeh = floor(rc_dt/dt);
switch rc_agg
    case 'first'
        timeVeh = 0;
    case 'middle'
        timeVeh = rc_dt/2;
    case 'last'
        timeVeh = rc_dt;
    case 'Null'
        for d_index=1:totDest
            arr_map_d(d_index);
        end
        TF=[];
        return;
    case 'inst'
        for d_index=1:totDest
            d=destinations(d_index);
            netCostMatrix=sparse(endN,strN,simTT(:,1),totNodes,totNodes);
            [parent, ~] = dijkstra(netCostMatrix, d);
            for n=1:totNodes
                incomingLinks = find(endN==n);
                outgoingLinks = find(strN==n);
                TF{n,1,d_index}=zeros(max(1,length(incomingLinks)),length(outgoingLinks));
                TF{n,1,d_index}(:,endN(outgoingLinks)==parent(n))=1;
            end
        end
        gap = [];
        return;
        
end
tVeh = floor(timeVeh/dt);

for d_index=1:totDest
    %find the arrival map and update the gap function
    [~,parent] = arr_map_d(d_index);
        
    %compute update of the turning fractions based on the projected
    %gradient
    for n=1:totNodes
        next_rc=1;
        incomingLinks = find(endN==n);
        outgoingLinks = find(strN==n);
        
        if length(outgoingLinks)<=1
            %no update is required for a non divergent node
            for t=1:totT
                TF{n,t,d_index}=ones(max(1,length(incomingLinks)),max(1,length(outgoingLinks)));
            end
        else
            for t=1:totT
                %update all turning fractions within the route choice
                %interval by same value
                if timeSteps(t)>=timeRC(min(length(timeRC),next_rc))
                    next_rc = next_rc+1;
                    %find the parent value
                    par = parent(n,min(totT+1,t+tVeh));
                    act_t(min(totT+1,t+tVeh))=true;
                end
                %updating of the turning fractions
                TF{n,t,d_index}=zeros(max(1,length(incomingLinks)),length(outgoingLinks));
                TF{n,t,d_index}(:,endN(outgoingLinks)==par)=1;
            end
        end
    end
    
    gap_dt=gap_dt+sum(sum(gap(:,2:end).*diff(cvn_up(:,:,d_index),1,2)));
    gap_rc=gap_rc+sum(sum(gap(:,act_t).*diff(cvn_up(:,[1,find(act_t)+gVeh-tVeh],d_index),1,2)));
end




%Nested function used for finding the destination based arrival map
    function [arr_map,parent] = arr_map_d(d_index)
        d=destinations(d_index);
        netCostMatrix=sparse(endN,strN,simTT(:,end),totNodes,totNodes);
        [par, dist] = dijkstra(netCostMatrix, d);
        parent = zeros(totNodes,totT+1);
        parent(:,totT+1) = par;
        arr_map = zeros(totNodes,totT+1);
        arr_map(:,totT+1)=dist+dt*totT;
        for t=totT+1:-1:1
            for n=1:totNodes
                if any(n==destinations)
                    if n~=d
                        arr_map(n,t)=inf;
                    else
                        arr_map(n,t)=(t-1)*dt;
                    end
                    continue;
                end
                outgoingLinks = find(strN==n);
                arr = inf;
                for l=outgoingLinks'
                    time=timeSteps(t)+simTT(l,t);
                    if time>=timeSteps(end)
                        val=time-dt*totT+arr_map(endN(l),end);
                    else
                        t1 = min(totT+1,max(t+1,1+floor(time/dt)));
                        t2 = min(totT+1,t1+1);
                        val = arr_map(endN(l),t1,:)+max(0,(1+time/dt-t1))*(arr_map(endN(l),t2,:)-arr_map(endN(l),t1,:));
                    end
                    if cvn_up(l,t,d_index)>0
                        gap(l,t) = gap(l,t) + val;
                    end
                    if val<=arr
                        arr=val;
                        parent(n,t) = endN(l);
                    end
                end
                for l=outgoingLinks'
                    if cvn_up(l,t,d_index)>0
                        gap(l,t) = gap(l,t) - arr;
                    end
                end
                arr_map(n,t)=arr;
            end
        end
    end
end
