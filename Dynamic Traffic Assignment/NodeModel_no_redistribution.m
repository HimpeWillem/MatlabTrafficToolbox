function [turnFlow] = NodeModel_no_redistribution(incomingLinks,outgoingLinks,sendingFlow,turningFractions,receivingFlow,caps,timeInterval)
%Node model
%
%SYNTAX
%   [flowIncomingLinks,flowOutgoingLinks] = NodeModel(incomingLinks,outgoingLinks,sendingFlow,turningFractions,receivingFlow,caps,timeInterval)
%
%DESCRIPTION
%   returns the transfer flow over a node
%
%INPUTS
%   incomingLinks: index of the incoming links
%   outgoingLinks: index of the outgoing links
%   sendingFlow: the sending flow
%   turningFractions: turning fractions at the node
%   receivingFlow: the receiving flow
%   caps: capacities
%   timeInterval: updating interval

    
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


intLinksToDo = length(incomingLinks);
dblReceivingFraction = zeros(size(outgoingLinks));
blnRecalculateSendingFlow = zeros(size(incomingLinks));
dblSendingCapFlows = zeros(length(incomingLinks),length(outgoingLinks));
dblCapFlows = zeros(size(outgoingLinks));
dblCapacities = zeros(length(incomingLinks),length(outgoingLinks));
blnLinkIsDone = zeros(size(incomingLinks));
dblTransferFlows = zeros(length(incomingLinks),length(outgoingLinks));

%determine overflow factor of outgoing constraints:
for j=1:length(outgoingLinks)
    dblDummy1 = 0;
    for i=1:length(incomingLinks)
        dblDummy1 = dblDummy1 + turningFractions(i,j)*sendingFlow(i); %holds the total sending flow:receiving link j
    end
    if dblDummy1 > 0
        dblReceivingFraction(j) = receivingFlow(j) / dblDummy1;  %holds the fraction of the receiving flow of link j over the total sending flow:link j
    else
        dblReceivingFraction(j) = inf;                                %if there%s no sending flow:receiving link j
    end
end

%check which incoming links may be supply constrained and which are
%definitely demand constrained:
for i = 1:length(incomingLinks)
    for j=1:length(outgoingLinks)
        if dblReceivingFraction(j) < 1
            if turningFractions(i,j)*sendingFlow(i) > 0 
                blnRecalculateSendingFlow(i) = 1;      %1 if the receiving flow is exceeded for any link j:which link i sends flow
                break
            end
        end
    end
    if blnRecalculateSendingFlow(i) == 1
        dblDummy1 = caps(incomingLinks(i))*timeInterval/sendingFlow(i);
        for j = 1:length(outgoingLinks)
            dblSendingCapFlows(i,j) = turningFractions(i,j)*sendingFlow(i)*dblDummy1; %capacity flow from link i:link j
            dblCapFlows(j) = dblCapFlows(j) + dblSendingCapFlows(i,j);                 %total flow at capacity:link j
        end
    else                                                                    %if no receiving flow is exceeded, link i can send all its flow
        for j = 1:length(outgoingLinks)
            dblCapacities(i,j) = turningFractions(i,j)*sendingFlow(i);
        end
        blnLinkIsDone(i) = 1;
        intLinksToDo = intLinksToDo - 1;
    end
end
 
%if there is at least one active supply constraint: calculate flows
MostRestrictiveOutgoingLink = 0;
while intLinksToDo > 0
    if MostRestrictiveOutgoingLink ~= length(outgoingLinks)
        MostRestrictiveOutgoingLink = MostRestrictiveOutgoingLink + 1; %needed:change MostRestrictiveOutgoingLink in case all RF%s become > 1
    else
        MostRestrictiveOutgoingLink = 1;
    end
    MostRestrictiveConstraint = 1;
    for j = 1:length(outgoingLinks)
        if dblCapFlows(j) > 0.000000000001 %prevents rounding off errors
            if receivingFlow(j) >= 0
                dblReceivingFraction(j) = receivingFlow(j) / dblCapFlows(j);
            else
                dblReceivingFraction(j) = 1;
            end
        else
            dblReceivingFraction(j) = 1;
        end
        if dblReceivingFraction(j) < MostRestrictiveConstraint                %the smallest RF is stored in dblConstraint
            MostRestrictiveOutgoingLink = j;
            MostRestrictiveConstraint = dblReceivingFraction(j);
        end
    end
    intLinksInFF = zeros(0,1);
    for i = 1:length(incomingLinks)
        if blnLinkIsDone(i) == 0                                  %if this SLi is not done yet
            if turningFractions(i,MostRestrictiveOutgoingLink) > 0          %...and sends flow:most restrictive RLj
                if sendingFlow(i) <= caps(incomingLinks(i))*timeInterval*MostRestrictiveConstraint    %if flow is not constrained
                    intLinksInFF(size(intLinksInFF,1)+1) = i;                   %all not constrained links are stored
                    for j = 1:length(outgoingLinks)
                        dblCapacities(i,j) = turningFractions(i,j) * sendingFlow(i);
                    end
                    blnLinkIsDone(i) = 1;
                    intLinksToDo = intLinksToDo - 1;
                end
            end
        end
    end
    if length(intLinksInFF) > 0           %at least 1 SL was not constrained
        for j = 1:length(outgoingLinks)
            for k = 1:length(intLinksInFF)
                dblCapFlows(j) = dblCapFlows(j) - dblSendingCapFlows(intLinksInFF(k), j);    %in the end run, we don%t consider flows from link i anymore
                receivingFlow(j) = receivingFlow(j) - turningFractions(intLinksInFF(k),j) * sendingFlow(intLinksInFF(k));
            end
        end
    else                                  %no links in FF, MostRestrictiveConstraint is applied on all SLi towards MostRestrictiveOutgoingLink
        for i = 1:length(incomingLinks)
            if blnLinkIsDone(i) == 0 
                if turningFractions(i,MostRestrictiveOutgoingLink) > 0 
                    dblCapFlows(MostRestrictiveOutgoingLink) = 0; 
                    receivingFlow(MostRestrictiveOutgoingLink) = 0;
                    for j = 1:length(outgoingLinks)
                        dblCapacities(i, j) = dblSendingCapFlows(i, j) * MostRestrictiveConstraint;
                    end
                    blnLinkIsDone(i) = 1;
                    intLinksToDo = intLinksToDo - 1;
                end
            end
        end
    end
end

%incoming flows:
dblMaxConstraint = zeros(size(incomingLinks));
for i = 1:length(incomingLinks)
    dblMaxConstraint(i) = 1;
    for j = 1:length(outgoingLinks)
        if turningFractions(i,j) * sendingFlow(i) > 0 
            dblMaxConstraint(i) = min(dblMaxConstraint(i), dblCapacities(i, j) / (turningFractions(i,j) * sendingFlow(i)));
            %for each sending link, the maximum constraint is determined
        end
    end
    flowIncomingLinks(i) = dblMaxConstraint(i) * sendingFlow(i);  %the flow of incoming links is calculated by the max constraint
end
  
%outgoing flows:
for i = 1:length(incomingLinks)
    for j = 1:length(outgoingLinks)
        dblTransferFlows(i, j) = turningFractions(i,j) * flowIncomingLinks(i);
    end
end
for j = 1:length(outgoingLinks)
    flowOutgoingLinks(j) = 0;
    for i = 1:length(incomingLinks)
        flowOutgoingLinks(j) = flowOutgoingLinks(j) + dblTransferFlows(i, j);
    end
end

turnFlow = repmat(flowIncomingLinks',1,length(outgoingLinks)).*turningFractions;