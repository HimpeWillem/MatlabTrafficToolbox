function [cvn_up,cvn_down] = LTM_SC(nodes,links,origins,destinations,ODmatrix,dt,totT,TF,varargin)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%link transmission model for network loading                              %
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
cvn_up = zeros(totLinks,totT+1);
cvn_down = zeros(totLinks,totT+1);

%local rename link properties (for shorter code)
fromNodes = links.fromNode;
toNodes = links.toNode;
freeSpeeds = links.freeSpeed;
capacities = links.capacity;
kJams = links.kJam;
lengths = links.length;
wSpeeds = capacities./(kJams-capacities./freeSpeeds);

normalNodes = setdiff(nodes.id,[origins,destinations]);

if nargin==8
    Lmod = 'PQ';
    Nmod = 'OC';
elseif nargin == 9
    Lmod = varargin{1};
    Nmod = 'OC';
else
    Lmod = varargin{1};
    Nmod = varargin{2};
end
    
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
                switch Lmod
                    case 'VQ'
                        RF(l_index) = calculateReceivingFlow_VQ(l,t);
                    case 'HQ'
                        RF(l_index) = calculateReceivingFlow_HQ(l,t);
                    case 'PQ'
                        RF(l_index) = calculateReceivingFlow_FQ(l,t);
                end
            end
            
            %compute transfer flows with the NODE MODEL
            switch Nmod
                case 'OC'
                    TransferFlow = NodeModel(nbIn,nbOut,SF,TF{n,t-1},RF,capacities(incomingLinks)*dt);
                case 'DP'
                    TransferFlow = NodeModel_demand_prop(nbIn,nbOut,SF,TF{n,t-1},RF,capacities(incomingLinks)*dt);
                case 'NR'
                    TransferFlow = NodeModel_no_redistribution(incomingLinks,outgoingLinks,SF,TF{n,t-1},RF,capacities,dt);
            end
            
            %update CVN values
            cvn_down(incomingLinks,t)=cvn_down(incomingLinks,t-1)+sum(TransferFlow,2);
            cvn_up(outgoingLinks,t)=cvn_up(outgoingLinks,t-1)+sum(TransferFlow,1)';
        end
    end
    
    %DESTINATION NODES<----------------------------------------------------------------------------------------------------------------------------
    %this nested function goes over all destination nodes
    loadDestinationNodes(t);
end

    %All nested function follow below:

    %Nested function for finding sending flows
    function SF = calculateSendingFlow(l,t)
        SF = capacities(l)*dt;
        time = timeSlices(t)-lengths(l)/freeSpeeds(l);
        val = findCVN(cvn_up(l,:),time,timeSlices,dt);
        SF = min(SF,val-cvn_down(l,t-1));
    end

    %Nested function for finding receiving flows for a vertical queue
    function RF = calculateReceivingFlow_VQ(l,t)
        RF = capacities(l)*dt;
    end

    %Nested function for finding receiving flows for a horizontal queue
    function RF = calculateReceivingFlow_HQ(l,t)
        RF = capacities(l)*dt;
        val = cvn_down(l,t-1)+kJams(l)*lengths(l);
        RF=min(RF,val-cvn_up(l,t-1));
    end

    %Nested function for finding receiving flows for a physical queue
    function RF = calculateReceivingFlow_FQ(l,t)
        RF = capacities(l)*dt;
        time = timeSlices(t)-lengths(l)/wSpeeds(l);
        val = findCVN(cvn_down(l,:),time,timeSlices,dt)+kJams(l)*lengths(l);
        RF = min(RF,val-cvn_up(l,t-1));
    end

    %Nested function used for finding CVN values inbetween time slices
    function val = findCVN(cvn,time,timeSlices,dt)
        if time<=timeSlices(1)
            val=0;
            return;
        elseif time>=timeSlices(end)
            val=cvn(end);
            return;
        else
            t1=ceil(time/dt);
            t2=t1+1;
            val = cvn(t1)+(time/dt-t1+1)*(cvn(t2)-cvn(t1));
        end
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
                SF = TF{o,t-1}.*sum(ODmatrix(o_index,:,t-1))*dt;
                cvn_up(l,t)=cvn_up(l,t-1) + SF;
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
                SF = capacities(l)*dt;
                SF = min(SF,findCVN(cvn_up(l,:),timeSlices(t)-lengths(l)/freeSpeeds(l),timeSlices,dt)-cvn_down(l,t-1));
                cvn_down(l,t)=cvn_down(l,t-1) + SF;
            end
        end 
    end
end