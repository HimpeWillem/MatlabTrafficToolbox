function [updateTF,gap_dt,gap_rc] = projGradTF(nodes,links,destinations,simTT,cvn_up,dt,totT,rc_dt,rc_agg,TF,alpha)
    
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

updateTF = num2cell(ones(size(nodes.id,1),totT,totDest));
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
        updateTF=[];
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
fracVeh = timeVeh/dt-tVeh;

for d_index=1:totDest
    %find the arrival map and update the gap function
    arr_map = arr_map_d(d_index);
    
    %no updating of the turning fractions is required
    if strcmp(rc_agg,'Null')
        continue;
    end
    
    %compute update of the turning fractions based on the projected
    %gradient
    for n=1:totNodes
        next_rc=1;
        incomingLinks = find(endN==n);
        outgoingLinks = find(strN==n);
        
        if length(outgoingLinks)<=1
            %no update is required for a non divergent node
            for t=1:totT
                updateTF{n,t,d_index}=zeros(max(1,length(incomingLinks)),max(1,length(outgoingLinks)));
            end
        else
            for t=1:totT
                %update all turning fractions within the route choice
                %interval by same value
                if timeSteps(t)>=timeRC(next_rc)
                    next_rc = next_rc+1;
                    act_t(min(totT+1,t+tVeh))=true;
                    
                    %compute arrival travel time at downstream end of all outgoing links
                    time=timeSteps(t)+timeVeh+(1-fracVeh)*(simTT(outgoingLinks,min(totT+1,t+tVeh)))+fracVeh*simTT(outgoingLinks,min(totT+1,t+tVeh+1));
                    t1 = min(totT+1,1+floor(time/dt));
                    vec1=sub2ind([totNodes,totT+1],endN(outgoingLinks),t1);
                    t2 = min(totT+1,t1+1);
                    vec2=sub2ind([totNodes,totT+1],endN(outgoingLinks),t2);
                    frac=(1+time/dt-t1);
                    %compute cost difference
                    costDiff = arr_map(vec1)+frac.*(arr_map(vec2)-arr_map(vec1))-((1-fracVeh)*arr_map(n,min(totT+1,t+tVeh))+fracVeh*arr_map(n,min(totT+1,t+tVeh+1)));
                    costDiff(time>=timeSteps(end)) = time(time>=timeSteps(end))-dt*totT+arr_map(endN(outgoingLinks(time>=timeSteps(end))),end)-((1-fracVeh)*(arr_map(n,min(totT+1,t+tVeh)))+fracVeh*arr_map(n,min(totT+1,t+tVeh+1)));
                    %next project the difference in the feasible probability space
                    P=zeros(1,length(outgoingLinks));
                    %only update turn probability if there is a positive turning fraction
                    %towards a longer route
                    if sum(TF{n,t,d_index}(1,costDiff>eps))>eps
                        %actual projection with scaling of the cost
                        %difference
                        P=-min(alpha*costDiff',TF{n,t,d_index}(1,:));
                        %quasi-reduced projection
                        [~,i]=min(costDiff);
                        P(i)=0;
                        P(i)=-sum(P);
                    end
                end
                %updating of the turning fractions
                updateTF{n,t,d_index}=zeros(max(1,length(incomingLinks)),length(outgoingLinks));
                updateTF{n,t,d_index}=repmat(P,max(1,length(incomingLinks)),1);
            end
        end
    end
    
    gap_dt=gap_dt+sum(sum(gap(:,2:end).*diff(cvn_up,1,2)));
    gap_rc=gap_rc+sum(sum(gap(:,act_t).*diff(cvn_up(:,[1,find(act_t)+gVeh-tVeh],d_index),1,2)));
end

%Nested function used for finding the destination based arrival map
    function arr_map = arr_map_d(d_index)
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
                    if val<arr
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
