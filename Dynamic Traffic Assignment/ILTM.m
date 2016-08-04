function [cvn_up,cvn_down,con_up,con_down] = ILTM(node_prop,links,origins,destinations,ODmatrix,dt,totT,TF...
    ,varargin)
%#codegen

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%intelligent link transmission assignment procedure                       %
%                                                                         %
%iterative updating of grid points                                        %
%intelligent activation/deactivations of grid points                      %
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

%make data link and data node work
totLinks=length(links.id);
totNodes=length(node_prop.nbIncomingLinks);
totDestinations=length(destinations);

%local rename link properties
strN = links.fromNode;
endN = links.toNode;
len = links.length;
cap = links.capacity;
kJm = links.kJam;
vind=links.vf_index;
vRt=links.vf_ratio;
wind=links.vw_index;
wRt=links.vw_ratio;

%local rename node properties
inL = node_prop.incomingLinksList;
outL = node_prop.outgoingLinksList;
nin = node_prop.nbIncomingLinks;
nout = node_prop.nbOutgoingLinks;
posIn = node_prop.positionFirstIn;
posOut = node_prop.positionFirstOut;
nTF = node_prop.nbTF;
inTF_index = node_prop.incomingLinksTFindex;
outTF_index = node_prop.outgoingLinksTFindex;
posTF = node_prop.positionFirstTF;

if nargin==8
    cvn_up = zeros(totLinks,totDestinations,totT+1);
    cvn_down = zeros(totLinks,totDestinations,totT+1);
    con_up = false(totLinks,totT+1);
    con_down = false(totLinks,totT+1);
    nodes2update = false(totNodes,totT+1);
    nodes2update(origins,:) = true;
    max_it_ILTM = 5000;
    marg_comp = false;
    precDNL = 10^-8;
elseif nargin==13
    cvn_up = varargin{1};
    cvn_down = varargin{2};
    con_up = varargin{3};
    con_down = varargin{4};
    nodes2update = varargin{5};
    max_it_ILTM = 5000;
    marg_comp = true;
    precDNL = 10^-8;
elseif nargin==15
    cvn_up = varargin{1};
    cvn_down = varargin{2};
    con_up = varargin{3};
    con_down = varargin{4};
    nodes2update = varargin{5};
    max_it_ILTM = varargin{6};
    precDNL = varargin{7};
    marg_comp = true;
else
    display('problem with number of inputs');
    return;
end

%max size of network
%(can be set dynamically as well need to check for speed)s
%(good practice to fix vector sizes on beforehand)
maxIncomingLinks = max(nin);
maxOutgoingLinks = max(nout);

%allocate memory to all local variables
RF_down_cvn_db = zeros(totLinks,1);
SF_up_cvn_db = zeros(totLinks,totDestinations);

incomingLinks = zeros(maxIncomingLinks,1);
outgoingLinks = zeros(maxOutgoingLinks,1);
tot_sendingFlow = zeros(maxIncomingLinks,1);
temp_capacities = zeros(maxIncomingLinks,1);
receivingFlow = zeros(maxOutgoingLinks,1);
sendingFlow = zeros(maxIncomingLinks,totDestinations);
temp_sendingFlow = zeros(maxIncomingLinks,totDestinations);
temp_receivingFlow = zeros(maxOutgoingLinks,totDestinations);
outgoingFlow = zeros(maxOutgoingLinks,totDestinations);

origins_t = union(origins,find(nin==0)');
destinations_t = union(destinations,find(nout==0)');

deltaChange = zeros(totNodes,1);
sortedNodes = zeros(totNodes,1);

turningFractions = zeros(maxIncomingLinks,maxOutgoingLinks);
turningFlows = zeros(maxIncomingLinks,maxOutgoingLinks);

%forward implicit scheme
%go sequentially over each time step
mean_it_iltm = 0;
max_it_iltm = 0;
totNodes_updates = 0;

%recalculate the cumulatives for each node at time step t until
%there are no more changes
for t=2:totT+1
    if ~any(nodes2update(:,t))
        continue;
    end
    
    %initialization
    %all helper functions need to be initialized at run time only once
    RF_down_cvn_bool = true(totLinks,1);
    SF_up_cvn_bool = true(totLinks,1);
    
    %ORIGIN NODES<----------------------------------------------------------------------------------------------------------------------------
    %update origin nodes if needed before other nodes
    if any(nodes2update(origins,t))
        for n_index=1:length(origins)
            n = origins(n_index);
            if nodes2update(n,t)
                totNodes_updates=totNodes_updates+1;
                outgoingLinks(1:nout(n)) = outL(posOut(n):(posOut(n)+nout(n)-1));
                for l_index=1:nout(n)
                    l=outgoingLinks(l_index);
                    temp_sendingFlow(1,:) = cvn_up(l,:,t-1) + TF(posTF(n)+l_index-1,:,t-1).*ODmatrix(n_index,:,t-1)*dt;
                    if sum(abs(temp_sendingFlow(1,:)-cvn_up(l,:,t)))>precDNL
                        %if an origin is updated the future grid-points
                        %also need checking because solution is CVN
                        nodes2update(n,min(totT+1,t+1)) = true;
                        cvn_up(l,:,t)=temp_sendingFlow(1,:);
                        
                        if sum(cvn_up(l,:,t))-sum(cvn_up(l,:,t-1)) < cap(l)*dt
                            con_up(outgoingLinks(l_index),t)=false;
                        else
                            con_up(outgoingLinks(l_index),t)=true;
                        end
                        
                        if vind(l)==-1
                            nodes2update(endN(l),t)=true;
                        else
                            nodes2update(endN(l),min(totT+1,t-vind(l)-1:t-vind(l)))=true;
                        end
                    end
                end
            end
        end
    end
    
    %first iteration nodes are sorted in no specific order
    n_i=1;
    for n=1:totNodes
        if nodes2update(n,t) && ~any(n == origins_t) && ~any(n == destinations_t)
            sortedNodes(n_i) = n;
            n_i=n_i+1;
        end
    end
    
    %amount of nodes that still need updating
    current_totNodes=n_i-1;
    
    % ACTUAL I-LTM ITERATIONS <-------------------------------------------------------------------------------------------------------------------------------
    %only on normal nodes
    %during iterations activate/deactivate nodes for updating based on
    %dynamics until no more nodes need updating
    it = 0;
    while it<max_it_ILTM && current_totNodes>0
        %as long as there are nodes to be updated do an iteration
        it = it + 1;
        
        %total amount of node updates for the ILTM procedure
        totNodes_updates=totNodes_updates + current_totNodes;
        
        %go over all ordered nodes that need updating
        for nIndex=1:current_totNodes
            n=sortedNodes(nIndex);
            
            incomingLinks(1:nin(n)) = inL(posIn(n):(posIn(n)+nin(n)-1));
            outgoingLinks(1:nout(n)) = outL(posOut(n):(posOut(n)+nout(n)-1));
            
            deltaChange(n)=0;
            
            %STANDARD NODES<-----------------------------------------------------------------------------------------------------------------------------------------------------------
            %most function calls will end up here
            
            %CALCULATE THE SENDING FLOW<----------------------------------------------------------------------------------------------------------------------------
            %this is the maximum number of vehicles comming from the incoming links that wants to travel over a node within this time interval
            %multiple solution are possible
            %so far only the simple triangular FD solution is implemented here
            %alternatives are a combined quadratic linear FD or
            %hypocritical waiting because of intersections
            for l_index=1:nin(n)
                l=incomingLinks(l_index);
                
                %initialize helper function if this is the first
                %time the link is evaluated
                if SF_up_cvn_bool(l)
                    SF_up_cvn_bool(l) = false;
                    SF_up_cvn_db(l,:) = cvn_up(l,:,(max(1,t+vind(l))))*(1-vRt(l)) - cvn_down(l,:,t-1);
                    if vind(l)<-1
                        SF_up_cvn_db(l,:) = SF_up_cvn_db(l,:)+vRt(l)*cvn_up(l,:,max(1,t+vind(l)+1));
                    end
                end
                
                sendingFlow(l_index,:)=SF_up_cvn_db(l,:);
                if vind(l)==-1
                    sendingFlow(l_index,:)=sendingFlow(l_index,:)+vRt(l)*cvn_up(l,:,t);
                end
                sendingFlow(l_index,:)=max(0,sendingFlow(l_index,:));
            end
            
            %CALCULATE RECEIVING FLOW<------------------------------------------------------------------------------------------------------------------------------
            %this is the maximum number of vehicles that can flow into the outgoing links
            %gridlock is still possible in the end solution
            %in some cases it is better to switch to vertical
            %queues then by setting the jam density to infinity
            for l_index=1:nout(n)
                l=outgoingLinks(l_index);
                
                %initialize helper function if this is the first
                %time the link is evaluated
                if RF_down_cvn_bool(l)
                    RF_down_cvn_bool(l) = false;
                    RF_down_cvn_db(l) = sum(cvn_down(l,:,max(1,t+wind(l)))).*(1-wRt(l)) - sum(cvn_up(l,:,t-1)) + kJm(l)*len(l);
                    if wind(l)<-1
                        RF_down_cvn_db(l)=RF_down_cvn_db(l)+wRt(l)*sum(cvn_down(l,:,max(1,t+wind(l)+1)));
                    end
                end
                
                receivingFlow(l_index)=RF_down_cvn_db(l);
                if wind(l)==-1
                    receivingFlow(l_index)=receivingFlow(l_index)+wRt(l).*sum(cvn_down(l,:,t));
                end
                receivingFlow(l_index) = min(cap(l)*dt,receivingFlow(l_index));
            end
            
            %CALCULATE TURNING FRACTIONS<---------------------------------------------------------------------------------------------------------------------------
            %split rates of all flow that runs over the node
            if nout(n)==1
                %Straight forward merge
                turningFractions(1:nin(n),1) = 1;
                outgoingFlow(1,:) = 1;
            else
                %more complex nodes with multiple outgoing links
                pl = posTF(n);
                turningFlows(:)=0;
                for tf_ind=1:nTF(n)
                    turningFlows(inTF_index(pl),outTF_index(pl))=turningFlows(inTF_index(pl),outTF_index(pl))+sum(sendingFlow(inTF_index(pl),:).*TF(pl,:,t-1),2);
                    pl=pl+1;
                end
                
                %derive turning fractions from the turningFlows
                for l_index=1:nin(n)
                    turningFractions(l_index,1:nout(n)) = turningFlows(l_index,1:nout(n))./sum(turningFlows(l_index,1:nout(n)));
                end
                turningFractions(isnan(turningFractions))=0;
            end
            
            %NODE MODEL<---------------------------------------------------------------------------------------------------------------------------------------------------------------------------
            %warning: the script for node model is initialized with static
            %fixed number of incoming or outgoing links
            %be carefull if number of incoming or outgoing links is
            %bigger then 7
            
            %prepocessing input
            for l_index=1:nin(n)
                tot_sendingFlow(l_index) = min(sum(sendingFlow(l_index,:)),cap(incomingLinks(l_index))*dt);
                temp_capacities(l_index)=cap(incomingLinks(l_index));
            end
            %capacity proportional
            %note that the node model can be further optimized to work faster
            %might use the dynamic version (it is a bit slower though)
            %or copy model here and make it destination based (not hard)
            [result_TurnFlow] = NodeModel(nin(n),nout(n),tot_sendingFlow(1:nin(n)),turningFractions(1:nin(n),1:nout(n)),receivingFlow(1:nout(n)),temp_capacities(1:nin(n)));
            result_sendingFlow=sum(result_TurnFlow,2);
            
            %can be that turningFractions need to be adjusted because
            %the sendingflow is less turning fractions might change a
            %bit only with small update steps (over multiple time
            %intervals linearization)
            
            %TRANSFER FLOWS && ACTIVATION OF NODES FOR NEXT ITERATION<--------------------------------------------------------------------------
            %does the node need updating?
            act_bool = false;
            temp_receivingFlow(1:nout(n),:) = 0;
            %first handle the incoming links
            tf_index=posTF(n)-1;
            for l_index=1:nin(n)
                tf_index=tf_index+1;
                %look for the active incoming links (links that send flow)
                if any(sendingFlow(l_index,:)>0)
                    %check congestion
                    if result_sendingFlow(l_index) < tot_sendingFlow(l_index) || result_sendingFlow(l_index)==temp_capacities(l_index)*dt
                        con_down(incomingLinks(l_index),t)=true;
                        %find the sending flow for each destination
                        %using FIFO rules distribution is proportional
                        temp_sendingFlow(l_index,:)=sendingFlow(l_index,:)/sum(sendingFlow(l_index,:))*result_sendingFlow(l_index);
                    else
                        con_down(incomingLinks(l_index),t)=false;
                        temp_sendingFlow(l_index,:)=sendingFlow(l_index,:);
                    end
                    
                    %find out if change to the node change is
                    %significant fot this incoming link
                    act_l_bool = sum(abs(cvn_down(incomingLinks(l_index),:,t)-(cvn_down(incomingLinks(l_index),:,t-1)+temp_sendingFlow(l_index,:))))>precDNL;
                    act_bool = act_bool|act_l_bool;
                    
                    %check if upstream moving change needs to be calculated
                    if act_l_bool
                        if wind(incomingLinks(l_index))==-1
                            nodes2update(strN(incomingLinks(l_index)),t)=true;
                            deltaChange(strN(incomingLinks(l_index)))=deltaChange(strN(incomingLinks(l_index)))+wRt(incomingLinks(l_index))*sum(abs(cvn_down(incomingLinks(l_index),:,t)-(cvn_down(incomingLinks(l_index),:,t-1)+temp_sendingFlow(l_index,:))));
                        else
                            nodes2update(strN(incomingLinks(l_index)),min(totT+1,t-wind(incomingLinks(l_index))-1:t-wind(incomingLinks(l_index)))) = true;
                        end
                    end
                    
                    active_d = sendingFlow(l_index,:)>0;
                    
                    %distribute them over the outgoing links
                    if nout(n)==1
                        temp_receivingFlow(1,active_d)=temp_receivingFlow(1,active_d)+temp_sendingFlow(l_index,active_d);
                    else
                        for l_index2=1:nout(n)
                            temp_receivingFlow(l_index2,active_d)=temp_receivingFlow(l_index2,active_d)+temp_sendingFlow(l_index,active_d).*TF(tf_index,active_d,t-1);
                            tf_index=tf_index+1;
                        end
                        tf_index=tf_index-1;
                    end
                else
                    for l_index2=2:nout(n)
                        tf_index=tf_index+1;
                    end
                    %no flow on this link
                    temp_sendingFlow(l_index,:) = 0;
                    con_down(incomingLinks(l_index),t)=false;
                end
            end
            
            %now check the outgoing links
            for l_index=1:nout(n)
                %find destination with significant change
                %then check for reactivation
                act_l_bool=sum(abs(cvn_up(outgoingLinks(l_index),:,t)-(cvn_up(outgoingLinks(l_index),:,t-1)+temp_receivingFlow(l_index,:))))>precDNL;
                act_bool = act_bool|act_l_bool;
                
                %check congestion
                if sum(temp_receivingFlow(l_index,:)) < receivingFlow(l_index)
                    con_up(outgoingLinks(l_index),t)=false;
                else
                    con_up(outgoingLinks(l_index),t)=true;
                end
                
                %check if downstream moving change needs to be calculated
                if act_l_bool
                    if vind(outgoingLinks(l_index))==-1
                        nodes2update(endN(outgoingLinks(l_index)),t)=true;
                        deltaChange(endN(outgoingLinks(l_index)))=deltaChange(endN(outgoingLinks(l_index)))+vRt(outgoingLinks(l_index))*sum(abs(cvn_up(outgoingLinks(l_index),:,t)-(cvn_up(outgoingLinks(l_index),:,t-1)+temp_receivingFlow(l_index,:))));
                    else
                        nodes2update(endN(outgoingLinks(l_index)),min(totT+1,t-vind(outgoingLinks(l_index))-1:t-vind(outgoingLinks(l_index))))=true;
                    end
                    
                    %check for changed flat CVN (instantenous) only for
                    %marginal updates needed
                    if marg_comp
                        pot_d=temp_receivingFlow(l_index,:)<=precDNL;
                        if any(pot_d) && sum(abs(cvn_up(outgoingLinks(l_index),pot_d,t)-(cvn_up(outgoingLinks(l_index),pot_d,t-1)+temp_receivingFlow(l_index,pot_d))))>precDNL
                            nodes2update(endN(outgoingLinks(l_index)),t:min(totT+1,t+1))=true;
                            deltaChange(endN(outgoingLinks(l_index)))=deltaChange(endN(outgoingLinks(l_index))) - dt*1/vind(outgoingLinks(l_index));
                        end
                    end
                    
                end
            end
            
            %only update all CVN if there is a significant change
            %on a node (nodes are always consistent because all CVN
            %are updated)
            if act_bool
                nodes2update(n,min(totT+1,t+1))=true;
            end
            for l_index=1:nin(n)
                cvn_down(incomingLinks(l_index),:,t)=cvn_down(incomingLinks(l_index),:,t-1)+temp_sendingFlow(l_index,:);
            end
            for l_index=1:nout(n)
                cvn_up(outgoingLinks(l_index),:,t)=cvn_up(outgoingLinks(l_index),:,t-1)+temp_receivingFlow(l_index,:);
            end
        end
        deltaChange(origins_t)=0;
        deltaChange(destinations_t)=0;
        %sort all changed nodes
        [~,sortedNodes(1:totNodes)] = sort(deltaChange(1:totNodes),1,'descend');
        %total active nodes for next iteration
        current_totNodes=sum(deltaChange(1:totNodes)>precDNL);
    end
    
    
    %DESTINATION NODES<----------------------------------------------------------------------------------------------------------------------------
    %are only updated if the upstream cumulative of
    %an incoming links is changed
    if any(nodes2update(destinations,t))
        for n_index=1:length(destinations)
            n = destinations(n_index);
            if nodes2update(n,t)
                totNodes_updates=totNodes_updates+1;
                incomingLinks(1:nin(n)) = inL(posIn(n):(posIn(n)+nin(n)-1));
                for l_index=1:nin(n)
                    l=incomingLinks(l_index);
                    RF_down_cvn_db(l)=inf;
                    %note that concidered intervals are always based on
                    %free flow travel time (destination is an infinit sink)
                    %IMPORTANT: this is a modelling assumption needs to be
                    %taken into account when setting up a network
                    %note that the cumulatives are found by looking
                    %into the up_cvn_ob of the previous iteration
                    %CALCULATE SENDING FLOW BASED ON FREE FLOW SPEED
                    %only n is needed this is the destination
                    %if multiple d are found an error should be thrown
                    temp_sendingFlow(1,:)=(1-vRt(l))*cvn_up(l,:,max(1,t+vind(l))) + vRt(l)*cvn_up(l,:,max(1,t+vind(l)+1));
                    if sum(abs(temp_sendingFlow(1,:)-cvn_down(l,:,t)))>precDNL
                        cvn_down(l,:,t)=temp_sendingFlow(1,:);
                    end
                end
            end
        end
    end
    mean_it_iltm = mean_it_iltm + max(1,it);
    max_it_iltm = max(max_it_iltm,it);
end

display(['average number of iterations: ',num2str(mean_it_iltm/totT)]);
display(['maximum number of iterations: ',num2str(max_it_iltm)]);
display(['total number of node updates: ',num2str(totNodes_updates)]);

end

