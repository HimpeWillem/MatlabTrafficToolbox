function [TF,gap_dt,gap_dt_s] = stochasticTF(nodes,links,destinations,simTT,cvn_up,dt,totT,rc_dt,rc_agg,theta)

totDest = length(destinations);
totNodes = length(nodes.id);
totLinks = length(links.id);
strN = links.fromNode;
endN = links.toNode;

if isempty(cvn_up)
    cvn_up=zeros(totLinks,totT+1,totDest);
end

updateTF = num2cell(ones(size(nodes.id,1),totT,totDest));
timeSteps = dt*[0:1:totT];
timeRC = rc_dt*[0:1:totT];
timeRC(timeRC>timeSteps(end))=[];

gap = zeros(totLinks,totT+1);
gap_dt = 0;
gap_s = zeros(totLinks,totT+1);
gap_dt_s = 0;
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
fracVeh = timeVeh/dt-tVeh;

for d_index=1:totDest
    %find the arrival map and update the gap function
    [arr_map,~] = arr_map_d(d_index);
       
    %use the arrival map order to compute the maximum perceived utility
    util_map = max_perc_util_d(arr_map,d_index);
    
    %compute update of the turning fractions based on the utiliy values
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
                    act_t(min(totT+1,t+tVeh))=true;

                    %compute the turn probabilities
                    P=zeros(1,length(outgoingLinks));
                    for i=1:length(outgoingLinks)
                        l=outgoingLinks(i);
                        if arr_map(n,min(totT+1,t+tVeh))>arr_map(endN(l),min(totT+1,t+tVeh))
                            time=timeSteps(min(totT+1,t+tVeh))+simTT(l,min(totT+1,t+tVeh));
                            n_down=endN(l);
                            if time>=timeSteps(end)
                                P(i)=exp((-simTT(l,min(totT+1,t+tVeh))+util_map(n_down,end)-util_map(n,min(totT+1,t+tVeh)))/theta);
                            else
                                t1 = min(totT+1,max(t+tVeh+1,1+floor(time/dt)));
                                t2 = min(totT+1,t1+1);
                                val = util_map(n_down,t1,:)+max(0,(1+time/dt-t1))*(util_map(n_down,t2,:)-util_map(n_down,t1,:));
                                P(i)=exp((-simTT(l,min(totT+1,t+tVeh))+val-util_map(n,min(totT+1,t+tVeh)))/theta);
                            end
                        end
                    end
                end
                %updating of the turning fractions
                TF{n,t,d_index}=zeros(max(1,length(incomingLinks)),length(outgoingLinks));
                TF{n,t,d_index}=repmat(P,max(1,length(incomingLinks)),1);
            end
        end
    end
    
    gap_dt=gap_dt+sum(sum(gap(:,2:end).*diff(cvn_up(:,:,d_index),1,2)));
    gap_dt_s=gap_dt_s+sum(sum(gap_s(:,2:end).*diff(cvn_up(:,:,d_index),1,2)));
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
                min_phi = inf;
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
                        phi = theta*log(cvn_up(l,t,d_index)-cvn_up(l,max(t-1,1),d_index))+val-(t-1)*dt;
                        gap_s(l,t) = phi;
                        if phi<=min_phi
                            min_phi=phi;
                        end
                    end
                    if val<=arr
                        arr=val;
                        parent(n,t) = endN(l);
                    end
                end
                for l=outgoingLinks'
                    if cvn_up(l,t,d_index)>0
                        gap(l,t) = gap(l,t) - arr;
                        gap_s(l,t) = gap_s(l,t)-min_phi;
                    end
                end
                arr_map(n,t)=arr;
            end
        end
    end

%Nested function used for finding the maximum perceived utility
    function [util_map] = max_perc_util_d(arr_map,d_index)
        d=destinations(d_index);
        util_map = zeros(totNodes,totT+1);
        %first do the last time slice
        [~,sorted_n]=sort(arr_map(:,end));
        for n_index=1:totNodes
            n=sorted_n(n_index);
            if any(n==destinations)
                if n~=d
                    util_map(n,t)=-inf;
                else
                    util_map(n,t)=0;
                end
                continue;
            end
            outgoingLinks = find(strN==n);
            for l=outgoingLinks'
                if arr_map(n,end)>arr_map(endN(l),end)
                    util_map(n,end)=util_map(n,end)+exp((-simTT(l,end)+util_map(endN(l),end))/theta);
                end
            end
            util_map(n,end) = theta*log(util_map(n,end));  
        end
        
        %next do the others in upwind order
        for t=totT:-1:1
            for n=1:totNodes
                if any(n==destinations)
                    if n~=d
                        util_map(n,t)=-inf;
                    else
                        util_map(n,t)=0;
                    end
                    continue;
                end
                outgoingLinks = find(strN==n);
                for l=outgoingLinks'
                    if arr_map(n,t)>arr_map(endN(l),t)
                        time=timeSteps(t)+simTT(l,t);
                        n_down=endN(l);
                        if time>=timeSteps(end)
                            util_map(n,t)=util_map(n,t)+exp((-simTT(l,t)+util_map(n_down,end))/theta);
                        else
                            t1 = min(totT+1,max(t+1,1+floor(time/dt)));
                            t2 = min(totT+1,t1+1);
                            val = util_map(n_down,t1,:)+max(0,(1+time/dt-t1))*(util_map(n_down,t2,:)-util_map(n_down,t1,:));
                            util_map(n,t)=util_map(n,t)+exp((-simTT(l,t)+val)/theta);
                        end
                    end
                end
                util_map(n,t) = theta*log(util_map(n,t));            
            end
        end
    end
end
