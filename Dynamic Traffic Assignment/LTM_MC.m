function [cvn_up,cvn_down] = LTM_MC(nodes,links,origins,destinations,ODmatrix,dt,totT,TF)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%link transmission assignment procedure                                   %
%                                                                         %
%destination based storing of commodities                                 %
%splitting rates at nodes based on TF                                     %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
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

%size of the network
totLinks = length(links.fromNode);
totDest = length(destinations);

%time slices for which a solution is build
timeSlices = [0:totT]*dt;

%cumulative vehicle numbers (cvn) are stored on both upstream and
%dowsntream link end of each link for every time slice
cvn_up = zeros(totLinks,totT+1,totDest);
cvn_down = zeros(totLinks,totT+1,totDest);

%local rename link properties (for shorter code)
fromNodes = links.fromNode;
toNodes = links.toNode;
freeSpeeds = links.freeSpeed;
capacities = links.capacity;
kJams = links.kJam;
lengths = links.length;
wSpeeds = capacities./(kJams-capacities./freeSpeeds);

normalNodes = setdiff(nodes.id,[origins,destinations]);

%forward explicit scheme
%go sequentially over each time step (first time step => all zeros)
for t=2:totT+1
    %ORIGIN NODES<--------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %this nested function goes over all origin nodes
    loadOriginNodes(t);
    
    %ACTUAL LTM <---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %go over all normal nodes in this time step
    for nIndex=1:length(normalNodes);
        %STANDARD NODES<--------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %most function calls will end up here
        n=normalNodes(nIndex);
                
        %CALCULATE THE SENDING FLOW<--------------------------------------------------------------------------------------------------------------------------------------------------------
        %this is the maximum number of vehicles comming from the
        %incoming links that want to travel over a node within this
        %time interval
        incomingLinks = find(toNodes==n);
        nbIn = length(incomingLinks);
        SF_d = zeros(nbIn,totDest);
        SF_tot = zeros(nbIn,1);
        SF = zeros(nbIn,1);
        for l_index=1:nbIn
            l=incomingLinks(l_index);
            SF_d(l_index,:) = calculateDestSendFlow(l,t);
            SF_tot(l_index) = sum(SF_d(l_index,:));
            SF(l_index) = min(capacities(l)*dt,SF_tot(l_index));
        end
        
        %CALCULATE RECEIVING FLOW<-----------------------------------------------------------------------------------------------------------------------------------------------------------
        %this is the maximum number of vehicles that can flow into the
        %outgoing links within this time interval
        outgoingLinks = find(fromNodes==n);
        nbOut = length(outgoingLinks);
        RF = zeros(nbOut,1);
        for l_index=1:nbOut
            l=outgoingLinks(l_index);
%             RF(l_index) = calculateReceivingFlow_VQ(l,t);
%             RF(l_index) = calculateReceivingFlow_HQ(l,t);
            RF(l_index) = calculateReceivingFlow_FQ(l,t);
        end
        
        %CALCULATE TURNING FRACTIONS<---------------------------------------------------------------------------------------------------------------------------
        %this calculates the split rates of all flows that run over the node
        TF_n = calculateTurningFractions(n,t-1);
         
        %compute transfer flows with the NODE MODEL
        TransferFlow = NodeModel(nbIn,nbOut,SF,TF_n,RF,capacities(incomingLinks)*dt);
         
        %update CVN values
        red = sum(TransferFlow,2)./(eps+SF_tot);
        for d = 1:totDest
            cvn_down(incomingLinks,t,d)=cvn_down(incomingLinks,t-1,d)+red.*SF_d(:,d);
            cvn_up(outgoingLinks,t,d)=cvn_up(outgoingLinks,t-1,d)+((red.*SF_d(:,d))'*TF{n,t-1,d})';
        end
    end
    
    %DESTINATION NODES<----------------------------------------------------------------------------------------------------------------------------
    %this nested function goes over all destination nodes
    loadDestinationNodes(t);
end

    %All nested function follow below:

    %Nested function for finding destination based sending flows
    function SF = calculateDestSendFlow(l,t)
        SFCAP = capacities(l)*dt;
        time = timeSlices(t)-lengths(l)/freeSpeeds(l);
        val = findCVN(cvn_up(l,:,:),time,timeSlices,dt);
        SF = val-cvn_down(l,t-1,:);
        if SF > SFCAP
            red = SFCAP/sum(SF);    
            SF = red*SF;
        end
        SF = max(0,SF);
    end

    %Nested function for finding receiving flows for a vertical queue
    function RF = calculateReceivingFlow_VQ(l,t)
        RF = capacities(l)*dt;
    end

    %Nested function for finding receiving flows for a horizontal queue
    function RF = calculateReceivingFlow_HQ(l,t)
        RF = capacities(l)*dt;
        val = sum(cvn_down(l,t-1,:),3)+kJams(l)*lengths(l)
        RF=min(RF,val-sum(cvn_up(l,t-1,:)));
    end

    %Nested function for finding receiving flows for a physical queue
    function RF = calculateReceivingFlow_FQ(l,t)
        RF = capacities(l)*dt;
        time = timeSlices(t)-lengths(l)/wSpeeds(l);
        val = findCVN(sum(cvn_down(l,:,:),3),time,timeSlices,dt)+kJams(l)*lengths(l);
        RF = min(RF,val-sum(cvn_up(l,t-1,:),3));
        RF = max(RF,0);
    end

    %Nested function for finding turning fractions
    function TF_n = calculateTurningFractions(n,t)
         %split rates of all flow that runs over the node
         if nbOut==1
             %Straight forward merge
             TF_n(1:nbIn,1) = 1;
         else
             %more complex nodes with multiple outgoing links
             TF_n=zeros(nbIn,nbOut);
             for d=1:totDest
                TF_n=TF_n+repmat(SF_d(:,d),1,nbOut).*TF{n,t,d};
             end
             %derive turning fractions from the turningFlows
             TF_n = TF_n./repmat(eps+sum(TF_n,2),1,nbOut);
         end
    end

    %Nested function that assigns the origin flow
    function loadOriginNodes(t)
        %update origin nodes
        for o_index=1:length(origins)
            o = origins(o_index);
            outgoingLinks = find(fromNodes==o);
            for l_index = 1:length(outgoingLinks)
                l=outgoingLinks(l_index);
                for d_index = 1:totDest
                    %calculation sending flow
                    SF_d = TF{o,t-1,d_index}.*sum(ODmatrix(o_index,d_index,t-1))*dt;
                    cvn_up(l,t,d_index)=cvn_up(l,t-1,d_index) + SF_d(l_index);
                end
            end
        end 
    end

    %Nested function that assigns the destination flow
    function loadDestinationNodes(t)
        %update origin nodes
        for d_index=1:length(destinations)
            d = destinations(d_index);
            incomingLinks = find(toNodes==d);
            for l_index=1:length(incomingLinks)
                l=incomingLinks(l_index);
                 for d_index = 1:totDest
                    %calculation sending flow
                    SF_d = findCVN(cvn_up(l,:,:),timeSlices(t)-lengths(l)/freeSpeeds(l),timeSlices,dt)-cvn_down(l,t-1,:);
                    cvn_down(l,t,:)=cvn_down(l,t-1,:) + SF_d;
                 end
            end
        end 
    end

    %Nested function used for finding CVN values inbetween time slices
    function val = findCVN(cvn,time,timeSlices,dt)
        if time<=timeSlices(1)
            val=cvn(1,1,:);
        elseif time>=timeSlices(end)
            val=cvn(1,end,:);
        else
            t1=ceil(time/dt);
            t2=t1+1;
            val = cvn(1,t1,:)+(time/dt-t1+1)*(cvn(1,t2,:)-cvn(1,t1,:));
        end
    end

end
