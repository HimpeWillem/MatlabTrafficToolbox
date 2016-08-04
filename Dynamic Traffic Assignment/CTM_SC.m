function [flows_up,flows_down,dens_n] = CTM_SC(nodes,links,origins,destinations,ODmatrix,dt,totT,TF)
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
totLinks=length(links.fromNode);

%time slices for which a solution is build
timeSlices = [0:totT]*dt;

%cumulative vehicle numbers (cvn) are stored on both upstream and
%dowsntream link end of each link for every time slice
dens_n = zeros(totLinks,totT+1);
flows_up = zeros(totLinks,totT);
flows_down = zeros(totLinks,totT);

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
    
    %ACTUAL CTM <---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
    %go over all normal nodes in this time step
    for nIndex=1:length(normalNodes);
        %STANDARD NODES<--------------------------------------------------------------------------------------------------------------------------------------------------------------------
        %most function calls will end up here
        n=normalNodes(nIndex);
        
        if isempty(TF{n,t-1}) || sum(sum((TF{n,t-1})))==0
            %no paths crossing this node
        else
            %active node
            
            %CALCULATE THE SENDING FLOW<--------------------------------------------------------------------------------------------------------------------------------------------------------
            %this is the maximum number of vehicles comming from the 
            %incoming links that want to travel over a node within this 
            %time interval
            incomingLinks = find(toNodes==n);
            nbIn = length(incomingLinks);
            SF = zeros(nbIn,1);
            for l_index=1:nbIn
                l=incomingLinks(l_index);
                SF(l_index) = calculateSendingFlow(l,t);
            end
            
            %CALCULATE RECEIVING FLOW<-----------------------------------------------------------------------------------------------------------------------------------------------------------
            %this is the maximum number of vehicles that can flow into the 
            %outgoing links within this time interval
            outgoingLinks = find(fromNodes==n);
            nbOut = length(outgoingLinks);
            RF = zeros(nbOut,1);
            for l_index=1:nbOut
                l=outgoingLinks(l_index);
%                 RF(l_index) = calculateReceivingFlow_VQ(l,t);
%                 RF(l_index) = calculateReceivingFlow_HQ(l,t);
                RF(l_index) = calculateReceivingFlow_FQ(l,t);
            end
            
            %compute transfer flows with the NODE MODEL
            TransferFlow = NodeModel(nbIn,nbOut,SF,TF{n,t-1},RF,capacities(incomingLinks));
%             TransferFlow = NodeModel_no_redistribution(incomingLinks,outgoingLinks,SF,TF{n,t-1},RF,capacities,dt);
            
            %update number of vehicles in a link
            flows_down(incomingLinks,t)=sum(TransferFlow,2);
            flows_up(outgoingLinks,t)=sum(TransferFlow,1)';
        end
    end
    
    %DESTINATION NODES<----------------------------------------------------------------------------------------------------------------------------
    %this nested function goes over all destination nodes
    loadDestinationNodes(t);
    
    
    %CTM update 
    dens_n(:,t) = dens_n(:,t-1)+1./lengths.*(flows_up(:,t)-flows_down(:,t))*dt;
end

    %All nested function follow below:

    %Nested function for finding sending flows
    function SF = calculateSendingFlow(l,t)
        SF = capacities(l);
        ki = dens_n(l,t-1);
        SF = min(SF,freeSpeeds(l)*ki);
    end

    %Nested function for finding receiving flows for a vertical queue
    function RF = calculateReceivingFlow_VQ(l,t)
        RF = capacities(l);
    end

    %Nested function for finding receiving flows for a horizontal queue
    function RF = calculateReceivingFlow_HQ(l,t)
        RF = capacities(l);
        RF=min(RF,wSpeeds(l)*kJams(l));
    end

    %Nested function for finding receiving flows for a physical queue
    function RF = calculateReceivingFlow_FQ(l,t)
        RF = capacities(l);
        ki = dens_n(l,t-1);
        RF = min(RF,wSpeeds(l)*(kJams(l)-ki));
    end

    %Nested function that assigns the origin flow
    function loadOriginNodes(t)
        %update origin nodes
        for o_index=1:length(origins)
            o = origins(o_index);
            outgoingLinks = find(fromNodes==o);
            for l_index=1:length(outgoingLinks)
                l=outgoingLinks(l_index);
                %calculation sending flow
                SF = TF{o,t-1}.*sum(ODmatrix(o_index,:,t-1));
                flows_up(l,t)=SF;
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
                %calculation sending flow
                SF = capacities(l);
                ki = dens_n(l,t-1);
                flows_down(l,t) = min(SF,freeSpeeds(l)*ki);
            end
        end 
    end
end